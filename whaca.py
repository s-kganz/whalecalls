# file io
from scipy.io import wavfile as wav
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
        self.NFFT = NFFT
        self.window = window

    def gen_spectro(self, file, process = True):
        # stub
        pass
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
    