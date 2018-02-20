
import numpy as np

def counts_error(vector):

    error_up = vector ** 0.5 
    error_down = vector ** 0.5 
 
    poisson_up=np.arange(31, dtype=np.float) 
    poisson_down=np.arange(31, dtype=np.float)

    poisson_up[0] = 1.841
    poisson_up[1] = 3.300
    poisson_up[2] = 4.638
    poisson_up[3] = 5.918
    poisson_up[4] = 7.163
    poisson_up[5] = 8.382
    poisson_up[6] = 9.584
    poisson_up[7] = 10.77
    poisson_up[8] = 11.95
    poisson_up[9] = 13.11
    poisson_up[10] = 14.27
    poisson_up[11] = 15.42
    poisson_up[12] = 16.56
    poisson_up[13] = 17.70
    poisson_up[14] = 18.83
    poisson_up[15] = 19.96
    poisson_up[16] = 21.08
    poisson_up[17] = 22.20
    poisson_up[18] = 23.32
    poisson_up[19] = 24.44
    poisson_up[20] = 25.55
    poisson_up[21] = 26.66
    poisson_up[22] = 27.76
    poisson_up[23] = 28.87
    poisson_up[24] = 29.97
    poisson_up[25] = 31.07
    poisson_up[26] = 32.16
    poisson_up[27] = 33.26
    poisson_up[28] = 34.35
    poisson_up[29] = 35.45
    poisson_up[30] = 36.54

    poisson_up = poisson_up - np.arange(31, dtype=np.float) 

    poisson_down[0] = 0.0
    poisson_down[1] = 0.173
    poisson_down[2] = 0.708
    poisson_down[3] = 1.367
    poisson_down[4] = 2.086
    poisson_down[5] = 2.840
    poisson_down[6] = 3.620
    poisson_down[7] = 4.419
    poisson_down[8] = 5.232
    poisson_down[9] = 6.057
    poisson_down[10] = 6.891
    poisson_down[11] = 7.734
    poisson_down[12] = 8.585
    poisson_down[13] = 9.441
    poisson_down[14] = 10.30
    poisson_down[15] = 11.17
    poisson_down[16] = 12.04
    poisson_down[17] = 12.92
    poisson_down[18] = 13.80
    poisson_down[19] = 14.68
    poisson_down[20] = 15.57
    poisson_down[21] = 16.45
    poisson_down[22] = 17.35
    poisson_down[23] = 18.24
    poisson_down[24] = 19.14
    poisson_down[25] = 20.03
    poisson_down[26] = 20.93
    poisson_down[27] = 21.84
    poisson_down[28] = 22.74
    poisson_down[29] = 23.65
    poisson_down[30] = 24.55

    poisson_down = np.arange(31, dtype=np.float) - poisson_down 

    error_up[(vector < 30.) & (vector > 0.)] = poisson_up[np.round(vector[(vector < 30.) & (vector > 0.)]).astype(int)] 
    error_down[(vector < 30.) & (vector > 0.)] = poisson_down[np.round(vector[(vector < 30.) & (vector > 0.)]).astype(int)] 

    error_up[vector < 0.] = poisson_up[0] 
    error_down[vector < 0.] = poisson_down[0]  

    return error_down, error_up 
