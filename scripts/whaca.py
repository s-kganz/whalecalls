# file io
from scipy.io import wavfile as wav
from obspy import read
# filter narrowband sound
from scipy.signal import medfilt
# specgram function
from matplotlib import mlab
import numpy as np

class whaca:

    DATA_TYPES = ["wav","mseed"]
    
    def __init__(self, db_thresh = 15, time_thresh = 0.5, width_thresh = 1000, NFFT = 512, window = None):
        # save threshold parameters
        self.db_thresh = db_thresh
        self.time_thresh = time_thresh
        self.width_thresh = width_thresh
        # save specgram parameters
        self.NFFT = NFFT
        self.window = window
    
    # file io
    
    def open_wav(self,f):
        self.rate,self.data = wav.read(f)
    
    def open_mseed(self,f):
        # TODO add support for multiple streams
        sound = read(f)[0]
        self.data = (sound.normalize() * (2 ** 31 - 1)).astype("int32")
        self.rate = sound.stats.sampling_rate
    
    # audio processing functions
    
    # subtract average intensity of each frequency
    def _sub_avg(self,s):
        new = s.copy()
        for r in new:
            r -= np.average(r)
        return new
    
    # eliminate broadband noises
    def _bb_reduc(self,s,freq):
        for i in range(0,np.ma.size(s,axis=1)):
            points = []
            for j, val in enumerate(s[:,i]):
                if val > self.db_thresh:
                    points.append(j)
            if len(points) != 0 and np.absolute(freq[max(points)] - freq[min(points)]) > self.width_thresh:
                s[:,i] = np.zeros(len(s[:,i])) + np.min(s)
        return s
    # generate audio specgram for call detection
    
    def gen_spectro(self, process = True):
        s,f,t = mlab.specgram(self.data,
                             NFFT = self.NFFT,
                             Fs = self.rate)
        # convert to db
        s = 10 * np.log10(s)
        s = np.flipud(s)
        if process:
            # TODO remove broadband sound
            s = self._bb_reduc(s,f)
            s = self._sub_avg(s)
            # restrict zeros
            s[s < 0] = 0
        return s,f,t
    
    # find sounds with given parameters
    
    def find_sounds(self):
        # stub
        # is this method necessary?
        pass
        
    # gets and sets
    
    def get_threshold_params(self):
        return {"db_thresh":self.db_thresh,
                "time_thresh":self.time_thresh,
                "width_thresh":self.width_thresh}
                
    def set_threshold_params(self, db_thresh = None, time_thresh = None, width_thresh = None):
        if db_thresh:
            self.db_thresh = db_thresh
        if time_thresh:
            self.time_thresh = time_thresh
        if width_thresh:
            self.width_thresh = width_thresh
            
    def set_spectro_params(self, NFFT = None, window = None):
        if NFFT:
            self.NFFT = NFFT
        if window:
            self.window = window
            
    def get_spectro_params(self):
        return {"NFFT":self.NFFT,
                "window":self.window}
    