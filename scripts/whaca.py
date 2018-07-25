# file io
from scipy.io import wavfile as wav
from obspy import read

# web resources
from urllib.request import urlretrieve

# filter narrowband sound
from scipy.signal import medfilt

# specgram function
from matplotlib import mlab
import numpy as np

import warnings


class whaca:
 
    def __init__(self, db_thresh=10, time_thresh=0.5, 
                 width_thresh=1000, NFFT=512, window=None):
        # save threshold parameters
        self.db_thresh = db_thresh
        self.time_thresh = time_thresh
        self.width_thresh = width_thresh
        # save specgram parameters
        self.NFFT = NFFT
        self.window = window
    
    # file io
    
    def open_wav(self, f):
        # call scipy module
        try:
            rate, data = wav.read(f)
        except ValueError:
            print("Cannot read wav. Maybe a bad url?")

        rate, data = wav.read(f)
        if data.ndim != 1:
            warnings.warn("Only mono files are supported, using left channel")
            data = data[:, 0]
        self.rate = rate
        self.data = data
            
    def open_mseed(self, f):
        # TODO add support for multiple streams
        try:
            sound = read(f)[0]
        except ValueError:
            print("Cannot read mseed. Maybe a bad url?")
            
        sound = read(f)[0]
        self.data = (sound.normalize() * (2 ** 31 - 1)).astype("int32")
        self.rate = sound.stats.sampling_rate
        
    def open_url(self, url, ext):
    
        data_types = {"wav": self.open_wav,
                      "mseed": self.open_mseed}
        # check if supported
        if ext in data_types.keys():
            f, headers = urlretrieve(url)
            data_types[ext](f)
        else:
            raise ValueError(ext + " is not a valid file type!")
    
    # audio processing functions
    
    # subtract average intensity of each frequency
    def _sub_avg(self, s):
        new = s.copy()
        for r in new:
            r -= np.average(r)
        return new
    
    # eliminate broadband noises
    def _bb_reduc(self, s, freq):
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
    
    # find longest sequence in array
    def _find_longest_sequence(self, arr):
        start = 0
        seq = [0, 0]
        for i in range(0, len(arr) - 1):
            if arr[i + 1] - arr[i] != 1 and i + 1 - start > len(seq):
                seq = arr[start:i+1]
                start = i + 1
        return seq[0], seq[-1]
      
    # generate audio specgram for call detection
    def gen_spectro(self, process=True):
        s, f, t = mlab.specgram(self.data,
                                NFFT=self.NFFT,
                                Fs=self.rate)
        # convert to db
        s = 10 * np.log10(s)
        s = np.flipud(s)
        if process:
            s = medfilt(s)
            s = self._bb_reduc(s, f)
            s = self._sub_avg(s)
            # restrict zeros
            s[s < 0] = 0
        return s, f, t
    
    # find sounds with given parameters
    
    def find_sounds(self):
        # stub
        # is this method necessary?
        pass
        
    # gets and sets
    
    def get_threshold_params(self):
        return {"db_thresh": self.db_thresh,
                "time_thresh": self.time_thresh,
                "width_thresh": self.width_thresh}
                
    def set_threshold_params(self, db_thresh=None, time_thresh=None, 
                             width_thresh=None):
        if db_thresh:
            self.db_thresh = db_thresh
        if time_thresh:
            self.time_thresh = time_thresh
        if width_thresh:
            self.width_thresh = width_thresh
            
    def set_spectro_params(self, NFFT=None, window=None):
        if NFFT:
            self.NFFT = NFFT
        if window:
            self.window = window
            
    def get_spectro_params(self):
        return {"NFFT": self.NFFT,
                "window": self.window}
    