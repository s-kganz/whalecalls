# file io
from scipy.io import wavfile as wav
from obspy import read

# web resources
from urllib.request import urlretrieve

# filter narrowband sound
from scipy.signal import medfilt

# call detection
from functools import reduce
import operator

# specgram function
from matplotlib import mlab
import numpy as np

import warnings


class whaca:
 
    def __init__(self, db_thresh=10, time_thresh=0.5,
                 width_thresh=1000, NFFT=512, window=None):
        '''
        Instantiate a new whaca class.

        Keyword arguments:
        db_thresh=10 : Intensity threshold for broadband noise reduction and
                       call filtering. (DB)
        time_thresh=0.5 : Time duration threshold for call filtering. (seconds)
        width_thresh=1000 : Bandwidth threshold for broadband noise reduction. (Hz)
        NFFT=512 : Sample size for Fourier transform on spectrogram generation.
                   Also roughly controls resolution of call detection.
        window=None : Window to be passed to Fourier transform function.

        '''
        # save threshold parameters
        self.db_thresh = db_thresh
        self.time_thresh = time_thresh
        self.width_thresh = width_thresh
        # save specgram parameters
        self.NFFT = NFFT
        self.window = window
    
    # file io
    
    def open_wav(self, f):
        '''Loads wav file data into class from the given file-like object or path, f.'''
        # call scipy module
        try:
            rate, data = wav.read(f)
        except ValueError:
            print("Cannot read wav. Maybe a bad url?")

        rate, data = wav.read(f)
        if data.ndim != 1:
            warnings.warn("Only mono files are supported, using left channel.")
            data = data[:, 0]
        self.rate = rate
        self.data = data
            
    def open_mseed(self, f):
        '''Loads mseed file data into class from the given file-like object or path, f.'''
        # TODO add support for multiple streams
        try:
            sound = read(f)[0]
        except ValueError:
            print("Cannot read mseed. Maybe a bad url?")
            
        sound = read(f)[0]
        self.data = (sound.normalize() * (2 ** 31 - 1)).astype("int32")
        self.rate = sound.stats.sampling_rate
        
    def open_url(self, url, ext):
        '''
        Loads data at specified url into class.

        Arguments:
        url: URL to be used.
        ext: Data type of the requested resource.
             Supported data types: 'wav', 'mseed'

        '''
        data_types = {"wav": self.open_wav,
                      "mseed": self.open_mseed}
        # check if supported
        if ext in data_types.keys():
            f, headers = urlretrieve(url)
            data_types[ext](f)
        else:
            raise ValueError(ext + " is not a valid file type!")
    
    # audio processing functions
    
    def _sub_avg(self, s):
        '''Subtracts average intensity of frequency bin from that bin in the given spectrogram.'''
        new = s.copy()
        for r in new:
            r -= np.average(r)
        return new
    
    def _bb_reduc(self, s, freq):
        '''
        Reduces broadband noise in the given spectrogram.

        Arguments:
        s: Spectrogram (2D numpy array) to be altered.
        freq: Frequency values corresponding to row of spectrogram.

        Returns:
        s: Altered spectrogram.

        '''
        # iterate along columns
        for i in range(0, np.ma.size(s, axis=1)):
            points = []
            for j, val in enumerate(s[:, i]):
                if val > self.db_thresh:
                    points.append(j)
            long = self._find_longest_sequence(points)
            if np.absolute(freq[long[1]] - freq[long[0]]) > self.width_thresh:
                s[:, i] = np.zeros(len(s[:, i])) + np.min(s)
        return s
    
    def _find_longest_sequence(self, arr):
        '''Finds longest consecutive sequence in array.'''
        start = 0
        seq = [0, 0]
        for i in range(0, len(arr) - 1):
            if arr[i + 1] - arr[i] != 1 and i + 1 - start > len(seq):
                seq = arr[start:i+1]
                start = i + 1
        return seq[0], seq[-1]
      
    def gen_spectro(self, process=True):
        s, f, t = mlab.specgram(self.data,
                                NFFT=self.NFFT,
                                Fs=self.rate)
        '''
        Generates spectrogram from data loaded into class.

        Keyword Arguments:
        process=True : Whether to apply median filter, broadband reduction,
                       and bin avg subtraction to spectrogram.

        Returns:
        s: Spectrogram produced from data.
        f: Frequency bins corresponding to rows of spectrogram.
        t: Timesteps corresponding to columns of spectrogram.
        '''
        # convert to db
        s = 10 * np.log10(s)
        s = np.flipud(s)
        if process:
            s = medfilt(s)
            s = self._bb_reduc(s, f)
            s = self._sub_avg(s)
            # restrict negative intensities
            s[s < 0] = 0
        return s, f, t
    
    def _group_consecutives(self, vals, step=1):
        '''Returns 2d list containing all consecutive subsets of array based off of given step size.'''
        run = []
        result = [run]
        expect = None
        for v in vals:
            if v == expect or expect is None:
                run.append(v)
            else:
                run = [v]
                result.append(run)
            expect = v + step
        return result
    
    def filter_sounds(self, spec, tstep):
        '''
        Filters spectrogram to only include potential whale calls.

        Arguments:
        spec: Spectrogram to be filtered.
        tstep: Difference in time between columns of spectrogram.

        Returns:
        spec: Filtered spectrogram.
        '''
        # drop values below threshold
        spec[spec < self.db_thresh] = 0
        # times containing valid intensities
        t = [col for col in range(0, np.ma.size(spec, axis=1))
             if any(i >= self.db_thresh for i in spec[:, col])]
        # keep time groups if they meet duration threshold
        groups = [lis for lis in self._group_consecutives(t)
                  if (len(lis) - 1) * tstep > self.time_thresh]
        # flatten list
        groups = reduce(operator.add, groups)
        # points not in groups need to be silenced
        silence = list(set(range(0, len(spec[0]))) - set(groups))
        spec[:, silence] = 0
        return spec
        
    # gets and sets
    
    def get_threshold_params(self):
        '''Returns dictionary of threshold parameters'''
        return {"db_thresh": self.db_thresh,
                "time_thresh": self.time_thresh,
                "width_thresh": self.width_thresh}
                
    def set_threshold_params(self, db_thresh=None, time_thresh=None,
                             width_thresh=None):
        '''Set threshold parameters db_thresh, time_thresh, and width_thresh.'''
        if db_thresh:
            self.db_thresh = db_thresh
        if time_thresh:
            self.time_thresh = time_thresh
        if width_thresh:
            self.width_thresh = width_thresh
            
    def set_spectro_params(self, NFFT=None, window=None):
        '''Set spectrogram parameters NFFT and window.'''
        if NFFT:
            self.NFFT = NFFT
        if window:
            self.window = window
            
    def get_spectro_params(self):
        '''Returns dictioary of threshold parmeters.'''
        return {"NFFT": self.NFFT,
                "window": self.window}
    