# -*- coding: utf-8 -*-
"""
Created on Mon Jun 13 16:45:08 2016

@author: manip.batm
"""

# -*- coding: utf-8 -*-
"""
Created on Wed Apr 13 15:15:32 2016

@author: manip.batm
"""

import re
import math
import numpy as np
import struct
import itertools
import pylab
import matplotlib.pyplot as plt
from operator import sub
from scipy import *
from jump import jump
from sweep_n_step2 import sweep_n_step2
from rabi import rabi
from multi_rabi_waitingtime import multi_rabi_waitingtime
from multi_rabi_amplitude import multi_rabi_amplitude


def step2_multiamplitude(X_A_start,Y_A_start,X_B_start,Y_B_start, X_C_start, Y_C_start,X_D_start,Y_D_start,X_E_start,Y_E_start,X_F_start,Y_F_start,X_G_start,Y_G_start,X_H_start,Y_H_start,X_I_start,Y_I_start,X_J_start,Y_J_start,X_A_stop,Y_A_stop,X_B_stop,Y_B_stop, X_C_stop, Y_C_stop,X_D_stop,Y_D_stop,X_E_stop,Y_E_stop,X_F_stop,Y_F_stop,X_G_stop,Y_G_stop,X_H_stop,Y_H_stop,X_I_stop,Y_I_stop,X_J_stop,Y_J_stop,wait_A,wait_AB,wait_B,wait_BC,wait_C,wait_CD,wait_D,wait_DE,wait_E,wait_EF,wait_F,wait_FG,wait_G,wait_GH,wait_H,wait_HI,wait_I,wait_IJ,wait_J,N_step,N_step2,delta_t, trigg_amp_init,trigg_amp_end, trigg_channel_length):
        
    Nindex=1
    ADC_samplingrate = 1.0e6;
    
#    AWG_samplingrate = 1.0e7;
       
    AWG_samplingrate = 1.0e9;
#    
    
    decimal=9    ####precision of the waiting time for example for a ns precision decimal=9 for 0.1us decimal=7
    trigg_length = 1* AWG_samplingrate/ADC_samplingrate;   # translate the sampling rate of the AWG in number of point related to the sampling rate of the ADC
    trigg_period = 1*trigg_length   # define the period of a trigg, related to the triggering of the ADC. With ADC trigg period=2*trigglength (one up and one down). With the LeCroy trigg period =1trigg length (just say to the LeCroy when it starts to measure, can do the same with the ADC)
    
   
    
    ###################################
    #########  active part
    ###################################
    

    wait_A_size = int(np.around(wait_A * ADC_samplingrate * trigg_period));
    wait_B_size = int(np.around(wait_B * ADC_samplingrate * trigg_period));
    wait_C_size = int(np.around(wait_C * ADC_samplingrate * trigg_period));
    wait_D_size = int(np.around(wait_D * ADC_samplingrate * trigg_period));
    wait_E_size = int(np.around(wait_E * ADC_samplingrate * trigg_period));
    wait_F_size = int(np.around(wait_F *  ADC_samplingrate * trigg_period));
    wait_G_size = int(np.around(wait_G * ADC_samplingrate * trigg_period));
    wait_H_size = int(np.around(wait_H * ADC_samplingrate * trigg_period));
    wait_I_size = int(np.around(wait_I * ADC_samplingrate * trigg_period));
    wait_J_size = int(np.around(wait_J * ADC_samplingrate * trigg_period));
    
    
    wait_AB_size = int(np.around(wait_AB * ADC_samplingrate * trigg_period));
    wait_BC_size = int(np.around(wait_BC * ADC_samplingrate * trigg_period));
    wait_CD_size = int(np.around(wait_CD * ADC_samplingrate * trigg_period));
    wait_DE_size = int(np.around(wait_DE * ADC_samplingrate * trigg_period));  
    wait_EF_size = int(np.around(wait_EF * ADC_samplingrate * trigg_period));
    wait_FG_size = int(np.around(wait_FG * ADC_samplingrate * trigg_period));   
    wait_GH_size = int(np.around(wait_GH * ADC_samplingrate * trigg_period));
    wait_HI_size = int(np.around(wait_HI * ADC_samplingrate * trigg_period));
    wait_IJ_size = int(np.around(wait_IJ * ADC_samplingrate * trigg_period));
        
        
    
     
    scan_length= int(wait_A_size+wait_B_size+wait_C_size+wait_D_size+wait_E_size+wait_F_size+wait_G_size+wait_H_size+wait_I_size+wait_J_size+wait_AB_size+wait_BC_size+wait_CD_size+wait_DE_size+wait_EF_size+wait_FG_size+wait_GH_size+wait_HI_size+wait_IJ_size) ;#length of a complete scan, with the ramp, steps, security ramp, etc
      
     
     
     
     
    ########################################
    #       INITIALISATION
    ########################################
    
    C1_DATA=np.zeros((Nindex+N_step*N_step2,scan_length));
    C2_DATA=np.zeros((Nindex+N_step*N_step2,scan_length));
    C3_DATA=np.zeros((Nindex+N_step*N_step2,scan_length));
    C4_DATA=np.zeros((Nindex+N_step*N_step2,scan_length));
    
    MARKER_C1_1=np.zeros((Nindex+N_step*N_step2,scan_length));
    MARKER_C2_1=np.zeros((Nindex+N_step*N_step2,scan_length));
    MARKER_C3_1=np.zeros((Nindex+N_step*N_step2,scan_length));
    MARKER_C4_1=np.ones((Nindex+N_step*N_step2,scan_length));
    
    MARKER_C1_2=np.zeros((Nindex+N_step*N_step2,scan_length));
    MARKER_C2_2=np.zeros((Nindex+N_step*N_step2,scan_length));
    MARKER_C3_2=np.zeros((Nindex+N_step*N_step2,scan_length));
    MARKER_C4_2=np.zeros((Nindex+N_step*N_step2,scan_length));
    
    JUMPS_DATA=np.zeros((Nindex+N_step*N_step2,4));
    SIZE=np.zeros((1,N_step*N_step2+Nindex))
    
    

    
    ######################################################################
    ###################         1=sweep          #########################
    ######################################################################
    #    
   
    ###########################################
    ############       X       ################
    ###########################################
    
    
    
    
    X_A_sweep=np.around(np.linspace(X_A_start,X_A_stop,N_step2),4)#
    
    X_B_sweep=np.around(np.linspace(X_B_start,X_B_stop,N_step2),4)#
    
    X_C_sweep=np.around(np.linspace(X_C_start,X_C_stop,N_step2),4)#
    
    X_D_sweep=np.around(np.linspace(X_D_start,X_D_stop,N_step2),4)#
    
    X_E_sweep=np.around(np.linspace(X_E_start,X_E_stop,N_step2),4)#
    
    X_F_sweep=np.around(np.linspace(X_F_start,X_F_stop,N_step2),4)#
    
    X_G_sweep=np.around(np.linspace(X_G_start,X_G_stop,N_step2),4)#
    
    X_H_sweep=np.around(np.linspace(X_H_start,X_H_stop,N_step2),4)#
    
    X_I_sweep=np.around(np.linspace(X_I_start,X_I_stop,N_step2),4)#
    
    X_J_sweep=np.around(np.linspace(X_J_start,X_J_stop,N_step2),4)#

    
    
    ###########################################
    ############       Y       ################
    ###########################################
    
    Y_A_sweep=np.around(np.linspace(Y_A_start,Y_A_stop,N_step2),4)#
    
    Y_B_sweep=np.around(np.linspace(Y_B_start,Y_B_stop,N_step2),4)#
    
    Y_C_sweep=np.around(np.linspace(Y_C_start,Y_C_stop,N_step2),4)#
    
    Y_D_sweep=np.around(np.linspace(Y_D_start,Y_D_stop,N_step2),4)#
    
    Y_E_sweep=np.around(np.linspace(Y_E_start,Y_E_stop,N_step2),4)#
    
    Y_F_sweep=np.around(np.linspace(Y_F_start,Y_F_stop,N_step2),4)#
    
    Y_G_sweep=np.around(np.linspace(Y_G_start,Y_G_stop,N_step2),4)#
    
    Y_H_sweep=np.around(np.linspace(Y_H_start,Y_H_stop,N_step2),4)#
    
    Y_I_sweep=np.around(np.linspace(Y_I_start,Y_I_stop,N_step2),4)#
    
    Y_J_sweep=np.around(np.linspace(Y_J_start,Y_J_stop,N_step2),4)#


      
    ###########################################
    ##########       trigg       ##############
    ###########################################  
      
    trigg_amp=np.around(np.linspace(trigg_amp_init,trigg_amp_end,N_step2),4) 
      
      
      
#    C1_DAT,C2_DAT,C3_DAT,C4_DAT,JUMPS_DAT,SIZE,M_C1_1,M_C2_1,M_C3_1,M_C4_1,M_C1_2,M_C2_2,M_C3_2,M_C4_2,AMP,OFFST=multi_rabi_amplitude(X_A_sweep[0],Y_A_sweep[0],X_B_sweep[0],Y_B_sweep[0], X_C_sweep[0],Y_C_sweep[0],X_D_sweep[0],Y_D_sweep[0],X_E_sweep[0],Y_E_sweep[0],X_F_sweep[0],Y_F_sweep[0],X_G_sweep[0],Y_G_sweep[0],X_H_sweep[0],Y_H_sweep[0],X_I_sweep[0],Y_I_sweep[0],X_J_sweep[0],Y_J_sweep[0],X_A_sweep[0],Y_A_sweep[0],X_B_sweep[0],Y_B_sweep[0], X_C_sweep[0],Y_C_sweep[0],X_D_sweep[0],Y_D_sweep[0],X_E_sweep[0],Y_E_sweep[0],X_F_sweep[0],Y_F_sweep[0],X_G_sweep[0],Y_G_sweep[0],X_H_sweep[0],Y_H_sweep[0],X_I_sweep[0],Y_I_sweep[0],X_J_sweep[0],Y_J_sweep[0],wait_A,wait_AB,wait_B,wait_BC, wait_C,wait_CD,wait_D,wait_DE,wait_E,wait_EF,wait_F,wait_FG,wait_G,wait_GH,wait_H,wait_HI,wait_I,wait_IJ,wait_J,wait_A,wait_AB,wait_B,wait_BC,wait_C,wait_CD,wait_D,wait_DE,wait_E,wait_EF,wait_F,wait_FG,wait_G,wait_GH,wait_H,wait_HI,wait_I,wait_IJ,wait_J,N_step,delta_t, trigg_amp_init, trigg_amp_end, trigg_channel_length)
    C1_DAT,C2_DAT,C3_DAT,C4_DAT,JUMPS_DAT,SIZE,M_C1_1,M_C2_1,M_C3_1,M_C4_1,M_C1_2,M_C2_2,M_C3_2,M_C4_2,AMP,OFFST=multi_rabi_amplitude(X_A_start,Y_A_start,X_B_start,Y_B_start, X_C_start, Y_C_start,X_D_start,Y_D_start,X_E_start,Y_E_start,X_F_start,Y_F_start,X_G_start,Y_G_start,X_H_start,Y_H_start,X_I_start,Y_I_start,X_J_start,Y_J_start,X_A_stop,Y_A_stop,X_B_stop,Y_B_stop, X_C_stop, Y_C_stop,X_D_stop,Y_D_stop,X_E_stop,Y_E_stop,X_F_stop,Y_F_stop,X_G_stop,Y_G_stop,X_H_stop,Y_H_stop,X_I_stop,Y_I_stop,X_J_stop,Y_J_stop,wait_A,wait_AB,wait_B,wait_BC, wait_C,wait_CD,wait_D,wait_DE,wait_E,wait_EF,wait_F,wait_FG,wait_G,wait_GH,wait_H,wait_HI,wait_I,wait_IJ,wait_J,wait_A,wait_AB,wait_B,wait_BC,wait_C,wait_CD,wait_D,wait_DE,wait_E,wait_EF,wait_F,wait_FG,wait_G,wait_GH,wait_H,wait_HI,wait_I,wait_IJ,wait_J,N_step,delta_t, trigg_amp[0], trigg_amp[0], trigg_channel_length)

    C1_DAT=np.delete(C1_DAT,(-1),axis=0)
    C2_DAT=np.delete(C2_DAT,(-1),axis=0)       
    C3_DAT=np.delete(C3_DAT,(-1),axis=0)       
    C4_DAT=np.delete(C4_DAT,(-1),axis=0)    
    
    M_C1_1=np.delete(M_C1_1,(-1),axis=0)       
    M_C2_1=np.delete(M_C2_1,(-1),axis=0)       
    M_C3_1=np.delete(M_C3_1,(-1),axis=0)       
    M_C4_1=np.delete(M_C4_1,(-1),axis=0)
    
    M_C1_2=np.delete(M_C1_2,(-1),axis=0)
    M_C2_2=np.delete(M_C2_2,(-1),axis=0)       
    M_C3_2=np.delete(M_C3_2,(-1),axis=0)       
    M_C4_2=np.delete(M_C4_2,(-1),axis=0) 
    
    
    
    
    C1_DATA[0:N_step+1][:]=C1_DAT
    C2_DATA[0:N_step+1][:]=C2_DAT
    C3_DATA[0:N_step+1][:]=C3_DAT
    C4_DATA[0:N_step+1][:]=C4_DAT    
    
    
    MARKER_C1_1[0:N_step+1][:]=M_C1_1
    MARKER_C2_1[0:N_step+1][:]=M_C2_1
    MARKER_C3_1[0:N_step+1][:]=M_C3_1
    MARKER_C4_1[0:N_step+1][:]=M_C4_1
      
    
    MARKER_C1_2[0:N_step+1][:]=M_C1_2
    MARKER_C2_2[0:N_step+1][:]=M_C2_2
    MARKER_C3_2[0:N_step+1][:]=M_C3_2
    MARKER_C4_2[0:N_step+1][:]=M_C4_2    
    
    
    C1_DAT=np.delete(C1_DAT,(0),axis=0)       
    C2_DAT=np.delete(C2_DAT,(0),axis=0)       
    C3_DAT=np.delete(C3_DAT,(0),axis=0)       
    C4_DAT=np.delete(C4_DAT,(0),axis=0)
    
    M_C1_1=np.delete(M_C1_1,(0),axis=0)      
    M_C2_1=np.delete(M_C2_1,(0),axis=0)       
    M_C3_1=np.delete(M_C3_1,(0),axis=0)       
    M_C4_1=np.delete(M_C4_1,(0),axis=0)
    
    M_C1_2=np.delete(M_C1_2,(0),axis=0)
    M_C2_2=np.delete(M_C2_2,(0),axis=0)       
    M_C3_2=np.delete(M_C3_2,(0),axis=0)       
    M_C4_2=np.delete(M_C4_2,(0),axis=0)
    
    
    for i in range(1,N_step2):
        
        C1_DAT,C2_DAT,C3_DAT,C4_DAT,JUMPS_DAT,SIZE,M_C1_1,M_C2_1,M_C3_1,M_C4_1,M_C1_2,M_C2_2,M_C3_2,M_C4_2,AMP,OFFST=multi_rabi_amplitude(X_A_start,Y_A_start,X_B_start,Y_B_start, X_C_start, Y_C_start,X_D_start,Y_D_start,X_E_start,Y_E_start,X_F_start,Y_F_start,X_G_start,Y_G_start,X_H_start,Y_H_start,X_I_start,Y_I_start,X_J_start,Y_J_start,X_A_stop,Y_A_stop,X_B_stop,Y_B_stop, X_C_stop, Y_C_stop,X_D_stop,Y_D_stop,X_E_stop,Y_E_stop,X_F_stop,Y_F_stop,X_G_stop,Y_G_stop,X_H_stop,Y_H_stop,X_I_stop,Y_I_stop,X_J_stop,Y_J_stop,wait_A,wait_AB,wait_B,wait_BC, wait_C,wait_CD,wait_D,wait_DE,wait_E,wait_EF,wait_F,wait_FG,wait_G,wait_GH,wait_H,wait_HI,wait_I,wait_IJ,wait_J,wait_A,wait_AB,wait_B,wait_BC,wait_C,wait_CD,wait_D,wait_DE,wait_E,wait_EF,wait_F,wait_FG,wait_G,wait_GH,wait_H,wait_HI,wait_I,wait_IJ,wait_J,N_step,delta_t, trigg_amp[i], trigg_amp[i], trigg_channel_length)
#        C1_DAT,C2_DAT,C3_DAT,C4_DAT,JUMPS_DAT,SIZE,M_C1_1,M_C2_1,M_C3_1,M_C4_1,M_C1_2,M_C2_2,M_C3_2,M_C4_2,AMP,OFFST=multi_rabi_amplitude(X_A_sweep[i],Y_A_sweep[i],X_B_sweep[i],Y_B_sweep[i], X_C_sweep[i],Y_C_sweep[i],X_D_sweep[i],Y_D_sweep[i],X_E_sweep[i],Y_E_sweep[i],X_F_sweep[i],Y_F_sweep[i],X_G_sweep[i],Y_G_sweep[i],X_H_sweep[i],Y_H_sweep[i],X_I_sweep[i],Y_I_sweep[i],X_J_sweep[i],Y_J_sweep[i],X_A_sweep[i],Y_A_sweep[i],X_B_sweep[i],Y_B_sweep[i], X_C_sweep[i],Y_C_sweep[i],X_D_sweep[i],Y_D_sweep[i],X_E_sweep[i],Y_E_sweep[i],X_F_sweep[i],Y_F_sweep[i],X_G_sweep[i],Y_G_sweep[i],X_H_sweep[i],Y_H_sweep[i],X_I_sweep[i],Y_I_sweep[i],X_J_sweep[i],Y_J_sweep[i],wait_A,wait_AB,wait_B,wait_BC, wait_C,wait_CD,wait_D,wait_DE,wait_E,wait_EF,wait_F,wait_FG,wait_G,wait_GH,wait_H,wait_HI,wait_I,wait_IJ,wait_J,wait_A,wait_AB,wait_B,wait_BC,wait_C,wait_CD,wait_D,wait_DE,wait_E,wait_EF,wait_F,wait_FG,wait_G,wait_GH,wait_H,wait_HI,wait_I,wait_IJ,wait_J,N_step,delta_t, trigg_amp_init, trigg_amp_end, trigg_channel_length)
    ##
        
        C1_DAT=np.delete(C1_DAT,(0),axis=0)
        C1_DAT=np.delete(C1_DAT,(-1),axis=0)
           
        C2_DAT=np.delete(C2_DAT,(0),axis=0)
        C2_DAT=np.delete(C2_DAT,(-1),axis=0)
           
        C3_DAT=np.delete(C3_DAT,(0),axis=0)
        C3_DAT=np.delete(C3_DAT,(-1),axis=0)
           
        C4_DAT=np.delete(C4_DAT,(0),axis=0)
        C4_DAT=np.delete(C4_DAT,(-1),axis=0)
        
    
        
        C1_DATA[1+i*N_step:N_step+i*N_step+1][:]=C1_DAT
        C2_DATA[1+i*N_step:N_step+i*N_step+1][:]=C2_DAT
        C3_DATA[1+i*N_step:N_step+i*N_step+1][:]=C3_DAT
        C4_DATA[1+i*N_step:N_step+i*N_step+1][:]=C4_DAT
        
    
    
        
        M_C1_1=np.delete(M_C1_1,(0),axis=0)
        M_C1_1=np.delete(M_C1_1,(-1),axis=0)
           
        M_C2_1=np.delete(M_C2_1,(0),axis=0)
        M_C2_1=np.delete(M_C2_1,(-1),axis=0)
           
        M_C3_1=np.delete(M_C3_1,(0),axis=0)
        M_C3_1=np.delete(M_C3_1,(-1),axis=0)
           
        M_C4_1=np.delete(M_C4_1,(0),axis=0)
        M_C4_1=np.delete(M_C4_1,(-1),axis=0)
        
        
        MARKER_C1_1[1+i*N_step:N_step+i*N_step+1][:]=M_C1_1
        MARKER_C2_1[1+i*N_step:N_step+i*N_step+1][:]=M_C2_1
        MARKER_C3_1[1+i*N_step:N_step+i*N_step+1][:]=M_C3_1
        MARKER_C4_1[1+i*N_step:N_step+i*N_step+1][:]=M_C4_1
        
    
        
        M_C1_2=np.delete(M_C1_2,(0),axis=0)
        M_C1_2=np.delete(M_C1_2,(-1),axis=0)
           
        M_C2_2=np.delete(M_C2_2,(0),axis=0)
        M_C2_2=np.delete(M_C2_2,(-1),axis=0)
           
        M_C3_2=np.delete(M_C3_2,(0),axis=0)
        M_C3_2=np.delete(M_C3_2,(-1),axis=0)
           
        M_C4_2=np.delete(M_C4_2,(0),axis=0)
        M_C4_2=np.delete(M_C4_2,(-1),axis=0)    
        
        
        MARKER_C1_2[1+i*N_step:N_step+i*N_step+1][:]=M_C1_2
        MARKER_C2_2[1+i*N_step:N_step+i*N_step+1][:]=M_C2_2
        MARKER_C3_2[1+i*N_step:N_step+i*N_step+1][:]=M_C3_2
        MARKER_C4_2[1+i*N_step:N_step+i*N_step+1][:]=M_C4_2
#        print 'ploppp de fin'
    
    
    
    
    
         
    SIZE[0][:]=np.shape(C4_DATA)[1];
    

    
     
    HIGHEST = array([C1_DATA.max(),C2_DATA.max(),C3_DATA.max(),C4_DATA.max()]);
    LOWEST  =array([C1_DATA.min(),C2_DATA.min(),C3_DATA.min(),C4_DATA.min()]);
    AMPLITUDE=np.zeros((1,4));
    for i in range(0,4):
        AMPLITUDE[0][i]=2*max(HIGHEST[i]-LOWEST[i],0.02)
    OFFSET = [0,0,0,0];

    # NORMALIZATION OF THE DATA TO BE TRANSMITTED TO THE AWG
#
    C1_DATA=2*(C1_DATA-OFFSET[0])/AMPLITUDE[0][0];
    C2_DATA=2*(C2_DATA-OFFSET[1])/AMPLITUDE[0][1];
    C3_DATA=2*(C3_DATA-OFFSET[2])/AMPLITUDE[0][2];
    C4_DATA=2*(C4_DATA-OFFSET[3])/AMPLITUDE[0][3];
#        
    return C1_DATA,C2_DATA,C3_DATA,C4_DATA,JUMPS_DATA,SIZE,MARKER_C1_1,MARKER_C2_1,MARKER_C3_1,MARKER_C4_1,MARKER_C1_2,MARKER_C2_2,MARKER_C3_2,MARKER_C4_2,AMPLITUDE,OFFSET

     