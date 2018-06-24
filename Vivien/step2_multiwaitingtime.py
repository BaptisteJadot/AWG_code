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


def step2_multiwaitingtime(X_A_start,Y_A_start,X_B_start,Y_B_start, X_C_start, Y_C_start,X_D_start,Y_D_start,X_E_start,Y_E_start,X_F_start,Y_F_start,X_G_start,Y_G_start,X_H_start,Y_H_start,X_I_start,Y_I_start,X_J_start,Y_J_start,X_A_stop,Y_A_stop,X_B_stop,Y_B_stop, X_C_stop, Y_C_stop,X_D_stop,Y_D_stop,X_E_stop,Y_E_stop,X_F_stop,Y_F_stop,X_G_stop,Y_G_stop,X_H_stop,Y_H_stop,X_I_stop,Y_I_stop,X_J_stop,Y_J_stop,wait_A_start,wait_AB_start,wait_B_start,wait_BC_start, wait_C_start,wait_CD_start,wait_D_start,wait_DE_start,wait_E_start,wait_EF_start,wait_F_start,wait_FG_start,wait_G_start,wait_GH_start,wait_H_start,wait_HI_start,wait_I_start,wait_IJ_start,wait_J_start,wait_A_stop,wait_AB_stop,wait_B_stop,wait_BC_stop,wait_C_stop,wait_CD_stop,wait_D_stop,wait_DE_stop,wait_E_stop,wait_EF_stop,wait_F_stop,wait_FG_stop,wait_G_stop,wait_GH_stop,wait_H_stop,wait_HI_stop,wait_I_stop,wait_IJ_stop,wait_J_stop,N_step,N_step2,delta_t_init,delta_t_fin, trigg_amp, trigg_channel_length):
    
    Nindex=1
    ADC_samplingrate = 1.0e6;
       
    AWG_samplingrate = 1.0e9;
    
    
    decimal=9    ####precision of the waiting time for example for a ns precision decimal=9 for 0.1us decimal=7
    trigg_length = 1* AWG_samplingrate/ADC_samplingrate;   # translate the sampling rate of the AWG in number of point related to the sampling rate of the ADC
    trigg_period = 1*trigg_length   # define the period of a trigg, related to the triggering of the ADC. With ADC trigg period=2*trigglength (one up and one down). With the LeCroy trigg period =1trigg length (just say to the LeCroy when it starts to measure, can do the same with the ADC)
    
    
    
    ###########################################
    #####  to test when change the amplitude
    ###########################################
    #db_atenuation_blue=db_atenuation_green=db_atenuation_red=141.7;
    #db_atenuation_black=157.4;
    #
    
    
    
#    db_atenuation_blue=db_atenuation_green=db_atenuation_red=89.1266;####6db
#
#    db_atenuation_black=100;####6db

    
    
     #X_A_start=X_H_start=-0.0*db_atenuation_blue/2;
#    X_B_start=-0.0001*db_atenuation_blue/2;
    #X_C_start=X_G_start=-0.006*db_atenuation_blue/2;
    #X_E_start=-0.006*db_atenuation_blue/2;
    #X_D_start=X_F_start=-0.023*db_atenuation_blue/2;
    #
    #Y_A_start=Y_H_start=-0.0212*db_atenuation_black/2;
#    Y_B_start=-0.0212*db_atenuation_black/2;
    #Y_C_start=Y_G_start=-0.0163*db_atenuation_black/2;
    #Y_E_start=-0.0135*db_atenuation_black/2;
    #Y_D_start=Y_F_start=-0.0015*db_atenuation_black/2;
    #
    #
    #
    #X_A_stop=X_H_stop=-0.0*db_atenuation_blue/2;
#    X_B_stop=-0.006*db_atenuation_blue/2;
    #X_C_stop=X_G_stop=-0.006*db_atenuation_blue/2;
    #X_E_stop=-0.006*db_atenuation_blue/2;
    #X_D_stop=X_F_stop=-0.023*db_atenuation_blue/2;
    #
    #Y_A_stop=Y_H_stop=-0.0212*db_atenuation_black/2;
#    Y_B_stop=-0.0212*db_atenuation_black/2;
    #Y_C_stop=Y_G_stop=-0.0163*db_atenuation_black/2;
    #Y_E_stop=-0.0163*db_atenuation_black/2;
    #Y_D_stop=Y_F_stop=-0.0015*db_atenuation_black/2;
    ##
    ###
    #N_step =51
    #N_step2 =25
    ###    
    #
    #
    #wait_A_start=12e-6
    #wait_AB_start=0e-9
    #wait_B_start=1e-9
    #wait_BC_start=500e-9
    #wait_C_start=2e-9
    #wait_CD_start=0e-9
    #wait_D_start=2e-9  #####
    #wait_DE_start=0e-9
    #wait_E_start=2e-9
    #wait_EF_start=500e-9
    #wait_F_start=1e-9
    #wait_FG_start=0e-9
    #wait_G_start=10e-9
    #wait_GH_start=0e-9
    #wait_H_start=10e-9
    #    
    #
    #
    #wait_A_stop=12e-6
    #wait_AB_stop=0e-9
    #wait_B_stop=1e-9
    #wait_BC_stop=500e-9
    #wait_C_stop=2e-9
    #wait_CD_stop=0e-9
    #wait_D_stop=52e-9  #####
    #wait_DE_stop=0e-9
    #wait_E_stop=2e-9
    #wait_EF_stop=500e-9
    #wait_F_stop=2e-9
    #wait_FG_stop=0e-9
    #wait_G_stop=10e-9
    #wait_GH_stop=0e-9
    #wait_H_stop=10e-9
    ###################################
    #########  active part
    ###################################
    
    deltat_A = (wait_A_stop-wait_A_start)
    deltat_AB = (wait_AB_stop-wait_AB_start)
    deltat_B = (wait_B_stop-wait_B_start)
    deltat_BC = (wait_BC_stop-wait_BC_start)
    deltat_C = (wait_C_stop-wait_C_start)
    deltat_CD = (wait_CD_stop-wait_CD_start)
    deltat_D = (wait_D_stop-wait_D_start)
    deltat_DE = (wait_DE_stop-wait_DE_start)
    deltat_E = (wait_E_stop-wait_E_start)
    deltat_EF = (wait_EF_stop-wait_EF_start)
    deltat_F = (wait_F_stop-wait_F_start)
    deltat_FG = (wait_FG_stop-wait_FG_start)
    deltat_G = (wait_G_stop-wait_G_start)
    deltat_GH = (wait_GH_stop-wait_GH_start)
    deltat_H = (wait_H_stop-wait_H_start)
    
    
    
    
    deltat=(wait_A_start+wait_B_start+wait_C_start+wait_D_start+wait_E_start+wait_AB_start+wait_BC_start+wait_CD_start+wait_DE_start+wait_EF_start+wait_F_start+wait_G_start+wait_FG_start+wait_H_start+wait_GH_start+wait_HI_start+wait_I_start+wait_IJ_start+wait_J_start)-(wait_A_stop+wait_B_stop+wait_C_stop+wait_D_stop+wait_E_stop+wait_AB_stop+wait_BC_stop+wait_CD_stop+wait_DE_stop+wait_F_stop+wait_EF_stop+wait_FG_stop+wait_G_stop+wait_GH_stop+wait_H_stop+wait_HI_stop+wait_I_stop+wait_IJ_stop+wait_J_stop)
    
    
    if deltat<0:    # =stop>start
        wait_A_size = int(np.around(wait_A_stop * ADC_samplingrate * trigg_period));
        wait_B_size = int(np.around(wait_B_stop * ADC_samplingrate * trigg_period));
        wait_C_size = int(np.around(wait_C_stop * ADC_samplingrate * trigg_period));
        wait_D_size = int(np.around(wait_D_stop * ADC_samplingrate * trigg_period));
        wait_E_size = int(np.around(wait_E_stop * ADC_samplingrate * trigg_period));
        wait_AB_size =int( np.around(wait_AB_stop * ADC_samplingrate * trigg_period));
        wait_BC_size =int( np.around(wait_BC_stop * ADC_samplingrate * trigg_period));
        wait_CD_size =int( np.around(wait_CD_stop * ADC_samplingrate * trigg_period));
        wait_DE_size =int( np.around(wait_DE_stop * ADC_samplingrate * trigg_period));
        wait_F_size = int(np.around(wait_F_stop * ADC_samplingrate * trigg_period));
        wait_EF_size =int( np.around(wait_EF_stop * ADC_samplingrate * trigg_period));
        wait_FG_size = int(np.around(wait_FG_stop * ADC_samplingrate * trigg_period));
        wait_G_size =int( np.around(wait_G_stop * ADC_samplingrate * trigg_period));
        wait_GH_size = int(np.around(wait_GH_stop * ADC_samplingrate * trigg_period));
        wait_H_size =int( np.around(wait_H_stop * ADC_samplingrate * trigg_period));
        wait_HI_size =int( np.around(wait_HI_start * ADC_samplingrate * trigg_period));
        wait_I_size =int( np.around(wait_I_start * ADC_samplingrate * trigg_period));
        wait_IJ_size = int(np.around(wait_IJ_stop * ADC_samplingrate * trigg_period));
        wait_J_size =int( np.around(wait_J_stop * ADC_samplingrate * trigg_period));
    else:
        wait_A_size = int(np.around(wait_A_start * ADC_samplingrate * trigg_period));
        wait_B_size = int(np.around(wait_B_start * ADC_samplingrate * trigg_period));
        wait_C_size = int(np.around(wait_C_start * ADC_samplingrate * trigg_period));
        wait_D_size = int(np.around(wait_D_start * ADC_samplingrate * trigg_period));
        wait_E_size = int(np.around(wait_E_start * ADC_samplingrate * trigg_period));
        wait_AB_size =int( np.around(wait_AB_start * ADC_samplingrate * trigg_period));
        wait_BC_size =int( np.around(wait_BC_start * ADC_samplingrate * trigg_period));
        wait_CD_size =int( np.around(wait_CD_start * ADC_samplingrate * trigg_period));
        wait_DE_size =int( np.around(wait_DE_start * ADC_samplingrate * trigg_period));
        wait_F_size =int( np.around(wait_F_start * ADC_samplingrate * trigg_period));
        wait_EF_size =int( np.around(wait_EF_start * ADC_samplingrate * trigg_period));
        wait_FG_size =int( np.around(wait_FG_start * ADC_samplingrate * trigg_period));
        wait_G_size =int( np.around(wait_G_start * ADC_samplingrate * trigg_period));
        wait_GH_size = int(np.around(wait_GH_stop * ADC_samplingrate * trigg_period));
        wait_H_size =int( np.around(wait_H_stop * ADC_samplingrate * trigg_period));
        wait_HI_size =int( np.around(wait_HI_start * ADC_samplingrate * trigg_period));
        wait_I_size =int( np.around(wait_I_start * ADC_samplingrate * trigg_period));
        wait_IJ_size = int(np.around(wait_IJ_stop * ADC_samplingrate * trigg_period));
        wait_J_size =int( np.around(wait_J_stop * ADC_samplingrate * trigg_period));
        
        
    
     
    scan_length= int(wait_A_size+wait_B_size+wait_C_size+wait_D_size+wait_E_size+wait_F_size+wait_G_size+wait_H_size+wait_I_size+wait_J_size+wait_AB_size+wait_BC_size+wait_CD_size+wait_DE_size+wait_EF_size+wait_FG_size+wait_GH_size+wait_HI_size+wait_IJ_size) ;#length of a complete scan, with the ramp, steps, security ramp, etc
      
     
     
     
     
    ########################################
    #       INITIALISATION
    ########################################
    #    print scan_length
    
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
    
    
    delta_A_X= (X_A_stop-X_A_start)
    delta_A_Y= (Y_A_stop-Y_A_start)
    
    delta_B_X= (X_B_stop-X_B_start)
    delta_B_Y= (Y_B_stop-Y_B_start)
    
    delta_C_X= (X_C_stop-X_C_start)
    delta_C_Y= (Y_C_stop-Y_C_start)
    
    delta_D_X= (X_D_stop-X_D_start)
    delta_D_Y= (Y_D_stop-Y_D_start)
    
    delta_E_X= (X_E_stop-X_E_start)
    delta_E_Y= (Y_E_stop-Y_E_start)
    
    delta_F_X= (X_F_stop-X_F_start)
    delta_F_Y= (Y_F_stop-Y_F_start)
    
    delta_G_X= (X_G_stop-X_G_start)
    delta_G_Y= (Y_G_stop-Y_G_start)
    
    delta_H_X= (X_H_stop-X_H_start)
    delta_H_Y= (Y_H_stop-Y_H_start)
    
    delta_I_X= (X_I_stop-X_I_start)
    delta_I_Y= (Y_I_stop-Y_I_start)
    
    delta_J_X= (X_J_stop-X_J_start)
    delta_J_Y= (Y_J_stop-Y_J_start)  
    
    
    delta_t_channel=delta_t_fin-delta_t_init
    
    #    
    #    
    #    ####  X  ####
    #    
    #  
    
      
    if delta_A_X==0:
       X_A=np.ones((N_step2))*X_A_start
    else:
        X_A=np.around(np.linspace(X_A_start,X_A_stop,N_step2),4)
    
    
    
    if delta_B_X==0:
       X_B=np.ones((N_step2))*X_B_start
    else:
        X_B=np.around(np.linspace(X_B_start,X_B_stop,N_step2),4)
    
    
    if delta_C_X==0:
       X_C=np.ones((N_step2))*X_C_start
    else:
        X_C=np.around(np.linspace(X_C_start,X_C_stop,N_step2),4)
    
    
    if delta_D_X==0:
       X_D=np.ones((N_step2))*X_D_start
    else:
        X_D=np.around(np.linspace(X_D_start,X_D_stop,N_step2),4)
        
    
    if delta_E_X==0:
       X_E=np.ones((N_step2))*X_E_start
    else:
        X_E=np.around(np.linspace(X_E_start,X_E_stop,N_step2),4)
        
        
    if delta_F_X==0:
       X_F=np.ones((N_step2))*X_F_start
    else:
        X_F=np.around(np.linspace(X_F_start,X_F_stop,N_step2),4)
        
    if delta_G_X==0:
       X_G=np.ones((N_step2))*X_G_start
    else:
        X_G=np.around(np.linspace(X_G_start,X_G_stop,N_step2),4)
    
    if delta_H_X==0:
       X_H=np.ones((N_step2))*X_H_start
    else:
        X_H=np.around(np.linspace(X_H_start,X_H_stop,N_step2),4)
  
  
    if delta_I_X==0:
       X_I=np.ones((N_step2))*X_I_start
    else:
        X_I=np.around(np.linspace(X_I_start,X_I_stop,N_step2),4)
  
    if delta_J_X==0:
       X_J=np.ones((N_step2))*X_J_start
    else:
        X_J=np.around(np.linspace(X_J_start,X_J_stop,N_step2),4)
    
    #    
    #    
    #    ####  Y  ####
    #    
    #  
    
    
    if delta_A_Y==0:
       Y_A=np.ones((N_step2))*Y_A_start
    else:
        Y_A=np.around(np.linspace(Y_A_start,Y_A_stop,N_step2),4)
    
    
    if delta_B_Y==0:
       Y_B=np.ones((N_step2))*Y_B_start
    else:
        Y_B=np.around(np.linspace(Y_B_start,Y_B_stop,N_step2),4)
    
    
    if delta_C_Y==0:
       Y_C=np.ones((N_step2))*Y_C_start
    else:
        Y_C=np.around(np.linspace(Y_C_start,Y_C_stop,N_step2),4)
    #    
    #    
    if delta_D_Y==0:
       Y_D=np.ones((N_step2))*Y_D_start
    else:
        Y_D=np.around(np.linspace(Y_D_start,Y_D_stop,N_step2),4)
    
    
    if delta_E_Y==0:
       Y_E=np.ones((N_step2))*Y_E_start
    else:
        Y_E=np.around(np.linspace(Y_E_start,Y_E_stop,N_step2),4)
        
        
    if delta_F_Y==0:
       Y_F=np.ones((N_step2))*Y_F_start
    else:
        Y_F=np.around(np.linspace(Y_F_start,Y_F_stop,N_step2),4)
        
    if delta_G_Y==0:
       Y_G=np.ones((N_step2))*Y_G_start
    else:
        Y_G=np.around(np.linspace(Y_G_start,Y_G_stop,N_step2),4)
    
    if delta_H_Y==0:
       Y_H=np.ones((N_step2))*Y_H_start
    else:
        Y_H=np.around(np.linspace(Y_H_start,Y_H_stop,N_step2),4)
    
    if delta_I_Y==0:
       Y_I=np.ones((N_step2))*Y_I_start
    else:
        Y_I=np.around(np.linspace(Y_I_start,Y_I_stop,N_step2),4)
          
    if delta_J_Y==0:
       Y_J=np.ones((N_step2))*Y_J_start
    else:
        Y_J=np.around(np.linspace(Y_J_start,Y_J_stop,N_step2),4)
        
        

    
    
    if delta_t_init-delta_t_fin == 0:
        delta_t=np.zeros((N_step2))
    else:
        delta_t=np.around(np.ones((N_step2))*np.arange(delta_t_init,delta_t_fin+delta_t_channel/N_step2,delta_t_channel/(N_step2-1)),decimal)

    
    ######################
    C1_DAT,C2_DAT,C3_DAT,C4_DAT,JUMPS_DAT,SIZE,M_C1_1,M_C2_1,M_C3_1,M_C4_1,M_C1_2,M_C2_2,M_C3_2,M_C4_2,AMP,OFFST=multi_rabi_waitingtime(X_A[0],Y_A[0],X_B[0],Y_B[0], X_C[0], Y_C[0],X_D[0],Y_D[0],X_E[0],Y_E[0],X_F[0],Y_F[0],X_G[0],Y_G[0],X_H[0],Y_H[0],X_I[0],Y_I[0],X_J[0],Y_J[0],wait_A_start,wait_AB_start,wait_B_start,wait_BC_start, wait_C_start,wait_CD_start,wait_D_start,wait_DE_start,wait_E_start,wait_EF_start,wait_F_start,wait_FG_start,wait_G_start,wait_GH_start,wait_H_start,wait_HI_start,wait_I_start,wait_IJ_start,wait_J_start,wait_A_stop,wait_AB_stop,wait_B_stop,wait_BC_stop,wait_C_stop,wait_CD_stop,wait_D_stop,wait_DE_stop,wait_E_stop,wait_EF_stop,wait_F_stop,wait_FG_stop,wait_G_stop,wait_GH_stop,wait_H_stop,wait_HI_stop,wait_I_stop,wait_IJ_stop,wait_J_stop,N_step,delta_t[0], trigg_amp, trigg_channel_length)
    
#    
  
    
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
        
#        print i        
    #
        C1_DAT,C2_DAT,C3_DAT,C4_DAT,JUMPS_DAT,SIZE,M_C1_1,M_C2_1,M_C3_1,M_C4_1,M_C1_2,M_C2_2,M_C3_2,M_C4_2,AMP,OFFST=multi_rabi_waitingtime(X_A[i],Y_A[i],X_B[i],Y_B[i], X_C[i], Y_C[i],X_D[i],Y_D[i],X_E[i],Y_E[i],X_F[i],Y_F[i],X_G[i],Y_G[i],X_H[i],Y_H[i],X_I[i],Y_I[i],X_J[i],Y_J[i],wait_A_start,wait_AB_start,wait_B_start,wait_BC_start, wait_C_start,wait_CD_start,wait_D_start,wait_DE_start,wait_E_start,wait_EF_start,wait_F_start,wait_FG_start,wait_G_start,wait_GH_start,wait_H_start,wait_HI_start,wait_I_start,wait_IJ_start,wait_J_start,wait_A_stop,wait_AB_stop,wait_B_stop,wait_BC_stop,wait_C_stop,wait_CD_stop,wait_D_stop,wait_DE_stop,wait_E_stop,wait_EF_stop,wait_F_stop,wait_FG_stop,wait_G_stop,wait_GH_stop,wait_H_stop,wait_HI_stop,wait_I_stop,wait_IJ_stop,wait_J_stop,N_step,delta_t[i], trigg_amp, trigg_channel_length)
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
        
        

#        
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
    
    
      
 
        
    return C1_DATA,C2_DATA,C3_DATA,C4_DATA,JUMPS_DATA,SIZE,MARKER_C1_1,MARKER_C2_1,MARKER_C3_1,MARKER_C4_1,MARKER_C1_2,MARKER_C2_2,MARKER_C3_2,MARKER_C4_2,AMPLITUDE,OFFSET
    
    
    
    
    
        
        
        
        
        
        
        
        
        
        
        
        
        
        
        
        