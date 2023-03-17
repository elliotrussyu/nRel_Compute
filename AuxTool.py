import time

# Self-defined rounding tool
def rnd(num :float, N :int)->float:
    """
    A simple rounding tool to round float num to a N decimals. 
    With >=5 round up and <= 4 round down.
    """
    mover = 10**N
    tempnum = int(num*mover)
    if tempnum%1 < 0.5:
        result = int(tempnum) / mover
    else:
        result = int(tempnum + 1) / mover
    return result


# A simple timer for performance testing
class Timer:
    def __init__(self):
        self.start_time = 0
        self.end_time = 0
        self.stamps = []
        self.tot_time = 0
        
    
    def reset(self):
        self.start_time = 0
        self.end_time = 0
        self.tot_time = 0
        self.stamps = []
        
    def begin(self):
        self.start_time = time.time()
    
    def stamp(self):
        self.stamps.append(time.time() - self.start_time)
        
    def finish(self):
        self.end_time = time.time()
        self.tot_time = rnd(self.end_time - self.start_time,6)
    
    def sh_runtime(self, pr = True):
        if pr: 
            print(self.tot_time)
        return self.tot_time
    
    def readable_time(self):
        pass
        # tot_time =  rndsd(end_time - start_time,6)
        # sec_decimals = tot_time % 1
    
    def avg_time(self, iterations):
        return rnd(self.tot_time / iterations , 6)