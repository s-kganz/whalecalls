# file io
from scipy.io import wavfile as wav
from obspy import read
import requests
# filter narrowband sound
from scipy.signal import medfilt
# specgram function
from matplotlib import mlab
import numpy as np

class whaca:
    
    def __init__(self, db_thresh = 10, time_thresh = 0.5, width_thresh = None, NFFT = 512, window = None):
        # save threshold parameters
        self.db_thresh = db_thresh
        self.time_thresh = time_thresh
        self.width_thresh = width_thresh
        # save specgram parameters
        self.NFFT = NFFT
        self.window = window
    
    def open_wav(self,f):
        # call scipy module
        self.rate,self.data = wav.read(f)
            
    def open_mseed(self,f):
        # TODO add support for multiple streams
        sound = read(f)[0]
        self.data = (sound.normalize() * (2 ** 31 - 1)).astype("int32")
        self.rate = sound.stats.sampling_rate
        
    def open_url(self,url):
    
        DATA_TYPES = {"wav":self.open_wav,
                  "mseed":self.open_mseed}
        
        # find last occurrence of dot
        di = url.rfind(".")
        # check the dot was found
        if di != -1:
            # determine file extension
            ext = url[di + 1:]
            # check if supported
            if ext in DATA_TYPES.keys():
                DATA_TYPES[ext](requests.get(url).content)
            else:
                raise ValueError(ext + " is not a valid file type!")
        else:
            raise ValueError(url + " does not look like a valid url!")   
            
    def gen_spectro(self, process = True):
        s,f,t = mlab.specgram(self.data,
                             NFFT = self.NFFT,
                             Fs = self.rate)
        # convert to db
        s = 10 * np.log10(s)
        s = np.flipud(s)
        return s,f,t
        
    def find_sounds(self):
        # stub
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
    