from __future__ import division
import numpy as np
from libs.SionAWG_class import SionAWG


def OneGateSpinMap(
        SweepChannel = 2,\
        AmplitudeStart = -4., \
        AmplitudeStop = -4.,\
        NumberPoints = 2,\
        PulseDuration = 100, \
        PulsePosition = 500,\
        WaveformDuration = 1000):

    """
    Returns NumberPoints pulses (PulseDuration,PulsePosition) on the gate
    SweepChannel with an amplitude varying from AmplitudeStart to
    AmplitudeStop.
    """

    norm = max(np.abs(AmplitudeStop),np.abs(AmplitudeStart))
    pulse_amplitudes = np.linspace(AmplitudeStart, AmplitudeStop, NumberPoints) / norm

    # Init the sequence
    sequence = {}
    sequence['WaitingSequenceElement'] = {}
    sequence['SequenceElements'] = {}
    sequence['Channels'] = {}
    sequence['NumberOfElements'] = NumberPoints

    # Build WaitingSequenceElement
    wait = {}
    wait['Name'] = 'Wait'
    wait['Size'] = WaveformDuration
    wait['Waveform'] = np.zeros((WaveformDuration,1))
    wait['Marker_1'] = []
    wait['Marker_2'] = []
    sequence['WaitingSequenceElement'] = wait

    # Build SequenceElements
    for i, pulse_amp in enumerate(pulse_amplitudes):
        i = i+1
        sequence['SequenceElements'][i] = {}
        sequence['SequenceElements'][i]['Index'] = i
        sequence['SequenceElements'][i]['Channels'] = {}
        for channel in range(1,5):
            if channel == SweepChannel:
                sequence['SequenceElements'][i]['Channels'][channel] = {}
                sequence['SequenceElements'][i]['Channels'][channel]['Name'] = 'C'+str(channel)+'P'+str(i)
                sequence['SequenceElements'][i]['Channels'][channel]['Size'] = WaveformDuration
                sequence['SequenceElements'][i]['Channels'][channel]['Waveform'] = np.zeros((WaveformDuration,1))
                sequence['SequenceElements'][i]['Channels'][channel]['Marker_1'] = []
                sequence['SequenceElements'][i]['Channels'][channel]['Marker_2'] = []
                sequence['SequenceElements'][i]['Channels'][channel]['Waveform'][PulsePosition:PulsePosition+PulseDuration]=pulse_amp
            else:
                sequence['SequenceElements'][i]['Channels'][channel] = wait

    # Set the Amplitude, Offset, Delay, Output of sequence['Channels']
    for channel in range(1,5):
        sequence['Channels'][channel] = {}
        sequence['Channels'][channel]['Offset'] = 0.0
        sequence['Channels'][channel]['Delay'] = 0.0
        sequence['Channels'][channel]['Output'] = False
        if channel == SweepChannel:
            amplitude = max(np.abs(AmplitudeStop), np.abs(AmplitudeStart))
            sequence['Channels'][channel]['Amplitude'] = max(amplitude,0.02)
        else:
            sequence['Channels'][channel]['Amplitude'] = 0.02

    return sequence
    
def SpinSeq_Twait(
        SweepChannel = 2,\
        AmplitudeStart = +2., \
        AmplitudeStop = +2.,\
        NumberPoints = 2,\
        PulseDuration = 100, \
        PulsePosition = 100,\
        WaveformDuration = 2000):

    """
    Returns NumberPoints pulses (PulseDuration,PulsePosition) on the gate
    SweepChannel with an amplitude varying from AmplitudeStart to
    AmplitudeStop.
    """

    norm = max(np.abs(AmplitudeStop),np.abs(AmplitudeStart))
    pulse_amplitudes = np.linspace(AmplitudeStart, AmplitudeStop, NumberPoints) / norm

    # Init the sequence
    sequence = {}
    sequence['WaitingSequenceElement'] = {}
    sequence['SequenceElements'] = {}
    sequence['Channels'] = {}
    sequence['NumberOfElements'] = NumberPoints

    # Build WaitingSequenceElement
    wait = {}
    wait['Name'] = 'Wait'
    wait['Size'] = WaveformDuration
    wait['Waveform'] = np.zeros((WaveformDuration,1))
    wait['Marker_1'] = []
    wait['Marker_2'] = []
    sequence['WaitingSequenceElement'] = wait

    # Build SequenceElements
    for i, pulse_amp in enumerate(pulse_amplitudes):
        i = i+1
        sequence['SequenceElements'][i] = {}
        sequence['SequenceElements'][i]['Index'] = i
        sequence['SequenceElements'][i]['Channels'] = {}
        for channel in range(1,5):
            if channel == SweepChannel:
                shape = np.concatenate((np.linspace(0.5,1.,501),np.linspace(1.,1.,101),np.linspace(1,0.5,501)))
                sequence['SequenceElements'][i]['Channels'][channel] = {}
                sequence['SequenceElements'][i]['Channels'][channel]['Name'] = 'C'+str(channel)+'P'+str(i)
                sequence['SequenceElements'][i]['Channels'][channel]['Size'] = WaveformDuration
                sequence['SequenceElements'][i]['Channels'][channel]['Waveform'] = np.zeros((WaveformDuration,1))
                sequence['SequenceElements'][i]['Channels'][channel]['Marker_1'] = []
                sequence['SequenceElements'][i]['Channels'][channel]['Marker_2'] = []
                sequence['SequenceElements'][i]['Channels'][channel]['Waveform'][PulsePosition:PulsePosition+len(shape),0]=shape
            else:
                sequence['SequenceElements'][i]['Channels'][channel] = wait

    # Set the Amplitude, Offset, Delay, Output of sequence['Channels']
    for channel in range(1,5):
        sequence['Channels'][channel] = {}
        sequence['Channels'][channel]['Offset'] = 0.0
        sequence['Channels'][channel]['Delay'] = 0.0
        sequence['Channels'][channel]['Output'] = False
        if channel == SweepChannel:
            amplitude = max(np.abs(AmplitudeStop), np.abs(AmplitudeStart))
            sequence['Channels'][channel]['Amplitude'] = max(amplitude,0.02)
        else:
            sequence['Channels'][channel]['Amplitude'] = 0.02

    return sequence

def OneGateSpinMapTwait(
        SweepChannel = 2,\
        Amplitude = -4, \
        PulseDurationStart = 5, \
        PulseDurationStop = 105, \
        NumberPoints = 20,\
        logscale = False,\
        WaveformDuration = 250):

    """
    Returns NumberPoints pulses (of amplitude Amplitude) on the gate
    SweepChannel with an duration varying from PulseDurationStart to
    PulseDurationStop. 
    """
    if logscale:
        pulse_durations = np.logspace(np.log10(PulseDurationStart), np.log10(PulseDurationStop), NumberPoints)
    else:
        pulse_durations = np.linspace(PulseDurationStart, PulseDurationStop, NumberPoints)

    # Init the sequence
    sequence = {}
    sequence['WaitingSequenceElement'] = {}
    sequence['SequenceElements'] = {}
    sequence['Channels'] = {}
    sequence['NumberOfElements'] = NumberPoints

    # Build WaitingSequenceElement
    wait = {}
    wait['Name'] = 'Wait'
    wait['Size'] = WaveformDuration
    wait['Waveform'] = np.zeros((WaveformDuration,1))
    wait['Marker_1'] = []
    wait['Marker_2'] = []
    sequence['WaitingSequenceElement'] = wait

    # Build SequenceElements
    for i, pulse_dur in enumerate(pulse_durations):
        i = i+1
        sequence['SequenceElements'][i] = {}
        sequence['SequenceElements'][i]['Index'] = i
        sequence['SequenceElements'][i]['Channels'] = {}
        for channel in range(1,5):
            if channel == SweepChannel:
                sequence['SequenceElements'][i]['Channels'][channel] = {}
                sequence['SequenceElements'][i]['Channels'][channel]['Name'] = 'C'+str(channel)+'T'+format(i,'02d')
                sequence['SequenceElements'][i]['Channels'][channel]['Size'] = WaveformDuration
                sequence['SequenceElements'][i]['Channels'][channel]['Waveform'] = np.zeros((WaveformDuration,1))
                sequence['SequenceElements'][i]['Channels'][channel]['Marker_1'] = []
                sequence['SequenceElements'][i]['Channels'][channel]['Marker_2'] = []
                sequence['SequenceElements'][i]['Channels'][channel]['Waveform'][-int(pulse_dur):] = -1.0
            else:
                sequence['SequenceElements'][i]['Channels'][channel] = wait

    # Set the Amplitude, Offset, Delay, Output of sequence['Channels']
    for channel in range(1,5):
        sequence['Channels'][channel] = {}
        sequence['Channels'][channel]['Offset'] = 0.0
        sequence['Channels'][channel]['Delay'] = 0.0
        sequence['Channels'][channel]['Output'] = False
        if channel == SweepChannel:
            amplitude = 1*abs(Amplitude)
            sequence['Channels'][channel]['Amplitude'] = max(amplitude,0.02)
        else:
            sequence['Channels'][channel]['Amplitude'] = 0.02

    return sequence
    
def PulseSequence():

    """
    Returns NumberPoints pulses (of amplitude Amplitude) on the gate
    SweepChannel with an duration varying from PulseDurationStart to
    PulseDurationStop. 
    """
    SweepChannel = 1
    Amplitudes = [+2.0,-2.0,0.0]
    TimingsStart = [1.0,1.0,1.0]
    TimingsStop = [1.0,1.0,1.0]
    NumberPoints = 5
    logscale = True
    WaveformDuration = 3
    SamplingRate = 10e6
    if len(Amplitudes) != len(TimingsStart) | len(Amplitudes) != len(TimingsStop):
        return 0
    
    Timings = np.zeros((NumberPoints,len(TimingsStart)))
    if logscale:
        for i,ti in enumerate(TimingsStart):
            Timings[:,i] = np.logspace(np.log10(TimingsStart[i]), np.log10(TimingsStop[i]), NumberPoints)
    else:
        for i,ti in enumerate(TimingsStart):
            Timings[:,i] = np.linspace(TimingsStart[i], TimingsStop[i], NumberPoints)
    
    # Init the sequence
    sequence = {}
    sequence['WaitingSequenceElement'] = {}
    sequence['SequenceElements'] = {}
    sequence['Channels'] = {}
    sequence['NumberOfElements'] = NumberPoints

    # Build WaitingSequenceElement
    WaitDuration = WaveformDuration
    wait = {}
    wait['Name'] = 'Wait'
    wait['Size'] = Time2Int(WaitDuration,SamplingRate)
    wait['Waveform'] = np.zeros((Time2Int(WaitDuration,SamplingRate),1))
    wait['Marker_1'] = []
    wait['Marker_2'] = []
    sequence['WaitingSequenceElement'] = wait

    # Build SequenceElements
    for i in range(1,NumberPoints+1):
        sequence['SequenceElements'][i] = {}
        sequence['SequenceElements'][i]['Index'] = i
        sequence['SequenceElements'][i]['Channels'] = {}
        for channel in range(1,5):
            if channel == SweepChannel:
                waveform = np.zeros((Time2Int(WaveformDuration,SamplingRate),1))
#                wfm = np.concatenate([[Ai*(i-1)*1.0/50.0 if Ai==-2 else Ai]*Time2Int(Timings[i-1,j],SamplingRate) for j,Ai in enumerate(Amplitudes)])
                wfm = np.concatenate([[Ai/2.0]*Time2Int(Timings[i-1,j],SamplingRate) for j,Ai in enumerate(Amplitudes)])
                waveform[-len(wfm):,0] = wfm
                sequence['SequenceElements'][i]['Channels'][channel] = {}
                sequence['SequenceElements'][i]['Channels'][channel]['Name'] = 'C'+str(channel)+'T'+format(i,'02d')
                sequence['SequenceElements'][i]['Channels'][channel]['Size'] = Time2Int(WaveformDuration,SamplingRate)
                sequence['SequenceElements'][i]['Channels'][channel]['Waveform'] = waveform
                sequence['SequenceElements'][i]['Channels'][channel]['Marker_1'] = []
                sequence['SequenceElements'][i]['Channels'][channel]['Marker_2'] = []
            else:
                sequence['SequenceElements'][i]['Channels'][channel] = wait

    # Set the Amplitude, Offset, Delay, Output of sequence['Channels']
    for channel in range(1,5):
        sequence['Channels'][channel] = {}
        sequence['Channels'][channel]['Offset'] = 0.0
        sequence['Channels'][channel]['Delay'] = 0.0
        sequence['Channels'][channel]['Output'] = False
        if channel == SweepChannel:
            amplitude = 2*max(Amplitudes)
            sequence['Channels'][channel]['Amplitude'] = max(amplitude,0.02)
        else:
            sequence['Channels'][channel]['Amplitude'] = 0.02

    return sequence
    
def Time2Int(t,SamplingRate):
    return int(t*1e-3*SamplingRate)
    
def BaptisteOnOff(
        SweepChannel = 3,\
        Amplitude = 0.5, \
        WaveformDuration = 1000):

    """
    Returns NumberPoints pulses (of amplitude Amplitude) on the gate
    SweepChannel with an duration varying from PulseDurationStart to
    PulseDurationStop. 
    """

    # Init the sequence
    sequence = {}
    sequence['WaitingSequenceElement'] = {}
    sequence['SequenceElements'] = {}
    sequence['Channels'] = {}
    sequence['NumberOfElements'] = 1

    # Build WaitingSequenceElement
    wait = {}
    wait['Name'] = 'Off'
    wait['Size'] = WaveformDuration
    wait['Waveform'] = np.zeros((WaveformDuration,1))
    wait['Marker_1'] = []
    wait['Marker_2'] = []
    sequence['WaitingSequenceElement'] = wait
    
    # Build SequenceElement
    for i in range(1,10):
        sequence['SequenceElements'][i] = {}
        sequence['SequenceElements'][i]['Index'] = 1
        sequence['SequenceElements'][i]['Channels'] = {}
        for channel in range(1,5):
            if channel == SweepChannel:
                sequence['SequenceElements'][i]['Channels'][channel] = {}
                sequence['SequenceElements'][i]['Channels'][channel]['Name'] = 'On'
                sequence['SequenceElements'][i]['Channels'][channel]['Size'] = WaveformDuration
                sequence['SequenceElements'][i]['Channels'][channel]['Waveform'] = -np.ones((WaveformDuration,1))
                sequence['SequenceElements'][i]['Channels'][channel]['Marker_1'] = []
                sequence['SequenceElements'][i]['Channels'][channel]['Marker_2'] = []
            else:
                sequence['SequenceElements'][i]['Channels'][channel] = wait
    
        # Set the Amplitude, Offset, Delay, Output of sequence['Channels']
        for channel in range(1,5):
            sequence['Channels'][channel] = {}
            sequence['Channels'][channel]['Offset'] = 0.0
            sequence['Channels'][channel]['Delay'] = 0.0
            sequence['Channels'][channel]['Output'] = False
            if channel == SweepChannel:
                amplitude = 2*Amplitude
                sequence['Channels'][channel]['Amplitude'] = max(amplitude,0.02)
            else:
                sequence['Channels'][channel]['Amplitude'] = 0.02
        return sequence
        
def PulsePosMap(
        SweepChannel = 4,\
        Amplitude = -2., \
        PulsePosStart = 1000, \
        PulsePosStop = 1000, \
        PulseDuration = 100, \
        NumberPoints = 1,\
        logscale = False,\
        WaveformDuration = 2000):

    """
    Returns NumberPoints pulses (of amplitude Amplitude) on the gate
    SweepChannel with an duration varying from PulseDurationStart to
    PulseDurationStop. 
    """
    pulse_pos = np.linspace(PulsePosStart, PulsePosStop, NumberPoints)

    # Init the sequence
    sequence = {}
    sequence['WaitingSequenceElement'] = {}
    sequence['SequenceElements'] = {}
    sequence['Channels'] = {}
    sequence['NumberOfElements'] = NumberPoints

    # Build WaitingSequenceElement
    wait = {}
    wait['Name'] = 'Wait'
    wait['Size'] = WaveformDuration
    wait['Waveform'] = np.zeros((WaveformDuration,1))
    wait['Marker_1'] = []
    wait['Marker_2'] = []
    sequence['WaitingSequenceElement'] = wait

    # Build SequenceElements
    for i, p0 in enumerate(pulse_pos):
        i = i+1
        sequence['SequenceElements'][i] = {}
        sequence['SequenceElements'][i]['Index'] = i
        sequence['SequenceElements'][i]['Channels'] = {}
        for channel in range(1,5):
            if channel == SweepChannel:
                sequence['SequenceElements'][i]['Channels'][channel] = {}
                sequence['SequenceElements'][i]['Channels'][channel]['Name'] = 'C'+str(channel)+'T'+format(i,'02d')
                sequence['SequenceElements'][i]['Channels'][channel]['Size'] = WaveformDuration
                sequence['SequenceElements'][i]['Channels'][channel]['Waveform'] = np.zeros((WaveformDuration,1))
                sequence['SequenceElements'][i]['Channels'][channel]['Marker_1'] = np.zeros((WaveformDuration,1))
                sequence['SequenceElements'][i]['Channels'][channel]['Marker_1'][1:120] = 1.    
                sequence['SequenceElements'][i]['Channels'][channel]['Marker_2'] = np.zeros((WaveformDuration,1))
                sequence['SequenceElements'][i]['Channels'][channel]['Marker_2'][1:120] = 1.
                sequence['SequenceElements'][i]['Channels'][channel]['Waveform'][p0:p0+PulseDuration] = 1.0
            else:
                sequence['SequenceElements'][i]['Channels'][channel] = wait

    # Set the Amplitude, Offset, Delay, Output of sequence['Channels']
    for channel in range(1,5):
        sequence['Channels'][channel] = {}
        sequence['Channels'][channel]['Offset'] = 0.0
        sequence['Channels'][channel]['Delay'] = 0.0
        sequence['Channels'][channel]['Output'] = False
        if channel == SweepChannel:
#            amplitude = 2*abs(Amplitude)
            amplitude = 2*abs(Amplitude)
            sequence['Channels'][channel]['Amplitude'] = max(amplitude,0.02)
        else:
            sequence['Channels'][channel]['Amplitude'] = 0.02

    return sequence
    
def MultiPulse(
    SweepChannel = 2,\
    Amplitude = 0.5, \
    PulseNumber = 10, \
    PulseDuration = 1000, \
    NumberPoints = 1,\
    WaveformDuration = 21000):

    """
    Returns NumberPoints pulses (of amplitude Amplitude) on the gate
    SweepChannel with an duration varying from PulseDurationStart to
    PulseDurationStop. 
    """
    ampl = np.linspace(1, 1, NumberPoints)

    # Init the sequence
    sequence = {}
    sequence['WaitingSequenceElement'] = {}
    sequence['SequenceElements'] = {}
    sequence['Channels'] = {}
    sequence['NumberOfElements'] = NumberPoints

    # Build WaitingSequenceElement
    wait = {}
    wait['Name'] = 'Wait'
    wait['Size'] = WaveformDuration
    wait['Waveform'] = np.zeros((WaveformDuration,1))
    wait['Marker_1'] = []
    wait['Marker_2'] = []
    sequence['WaitingSequenceElement'] = wait

    wfm_shape = np.array(([0]*PulseDuration+[-1]*PulseDuration)*PulseNumber+[0]*(WaveformDuration-PulseNumber*2*PulseDuration))
    wfm_shape = wfm_shape.reshape((WaveformDuration,1))

    # Build SequenceElements
    for i, ai in enumerate(ampl):
        i = i+1
        sequence['SequenceElements'][i] = {}
        sequence['SequenceElements'][i]['Index'] = i
        sequence['SequenceElements'][i]['Channels'] = {}
        for channel in range(1,5):
            if channel == SweepChannel:
                sequence['SequenceElements'][i]['Channels'][channel] = {}
                sequence['SequenceElements'][i]['Channels'][channel]['Name'] = 'C'+str(channel)+'A'+format(i,'02d')
                sequence['SequenceElements'][i]['Channels'][channel]['Size'] = WaveformDuration
                sequence['SequenceElements'][i]['Channels'][channel]['Waveform'] = wfm_shape*ai
                sequence['SequenceElements'][i]['Channels'][channel]['Marker_1'] = np.zeros((WaveformDuration,1))
                sequence['SequenceElements'][i]['Channels'][channel]['Marker_1'][1:100] = 1.    
                sequence['SequenceElements'][i]['Channels'][channel]['Marker_2'] = np.zeros((WaveformDuration,1))
                sequence['SequenceElements'][i]['Channels'][channel]['Marker_2'][1:100] = 1.
            else:
                sequence['SequenceElements'][i]['Channels'][channel] = wait

    # Set the Amplitude, Offset, Delay, Output of sequence['Channels']
    for channel in range(1,5):
        sequence['Channels'][channel] = {}
        sequence['Channels'][channel]['Offset'] = 0.0
        sequence['Channels'][channel]['Delay'] = 0.0
        sequence['Channels'][channel]['Output'] = False
        if channel == SweepChannel:
#            amplitude = 2*abs(Amplitude)
            amplitude = 2*abs(Amplitude)
            sequence['Channels'][channel]['Amplitude'] = max(amplitude,0.02)
        else:
            sequence['Channels'][channel]['Amplitude'] = 0.02

    return sequence
    
    
if __name__ == '__main__':
    sion = SionAWG('192.168.137.2', 4000)
#    seq = MultiPulse()
#    seq = BaptisteOnOff()
#    seq = SpinSeq_Twait()
    seq = OneGateSpinMap()
    sion.openCom()
    sion.DeleteAllWaveforms()
#    sion.SendSequenceLight(sequence = seq)
    sion.SendSequenceBaptiste(sequence = seq)
    sion.closeCom()
