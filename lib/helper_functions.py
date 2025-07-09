import numpy as np
import math
csc_veto = 0.978
dt_veto = 0.516
def weight_calc(llp_ct, new_ctau, old_ctau, nLLP = 2, flag = False):
    source = np.exp(-1.0*llp_ct/old_ctau)/old_ctau**nLLP
    weight = 1.0/new_ctau**nLLP * np.exp(-1.0*llp_ct/new_ctau)/source
    return weight


def signal_systematics(k, category, bins = 4):
    
    if "csccsc" in category:
        if "low" in k: 
            JES = [0.0189,0.0266,0.0240,0.0172]
            pileup = [0.0541, 0.0096, 0.0208, 0.0377]
            higgsPt = [ 0.0314,  0.0187,  0.0125, 0.0283, 0.0254,0.0151,0.0099,0.0232]
        elif "high" in k: 
            JES =  [0.0073,0.0043,0.0074,0.0183]
            pileup = [0.0112,0.0041,0.0516, 0.0171]
            higgsPt = [0.0111,0.0876,0.0795,0.0105,0.0092,0.0697,0.0632,0.0087]
        else: assert(False)
        dnn = [0.24]*bins
        timespread = [0.14]*bins
    elif "dtcsc" in category or "cscdt" in category:
        if "low" in k: 
            JES = [0.0076,0.0149,0.0159,0.0088]
            pileup = [0.01616,0.0390,0.0171,0.0164]
            higgsPt = [0.0343,0.0132,0.0125,0.0337,0.0285,0.0109,0.0101,0.0278]
        elif "high" in k: 
            JES =  [0.0216,0.0071,0.0046,0.0266]
            pileup = [0.0087,0.0131,0.0090,0.0067]
            higgsPt = [0.0102,0.0894,0.0875,0.0065,0.0084,0.0711,0.0698,0.0053]
        dnn = [0.12]*bins
        timespread = [0.07]*bins
        rpcBX = [0.077]*bins
        rpcMatch = [0.025]*bins
    if "4B" in k and "low" in k:
        csc_nhits = [0.15277692806274012, 0.14861595313954679, 0.19340720026134428, 0.1731532313851316]
    if "4B" in k and "high" in k:
        csc_nhits = [ 0.16458357643149102, 0.14429296134805625, 0.2272789539532527, 0.22895718300273327]

    if "4Tau" in k and "low" in k:
        csc_nhits = [0.22594125657443997, 0.153069776723737, 0.06996316860029084, 0.11404482922738701]

    if "4Tau" in k and "high" in k:
        csc_nhits = [ 0.2384315323776277, 0.2020066992989643, 0.1647089377667703, 0.1658265421718068]

    lumi = [0.02]*bins
    xsec = [0.067]*bins+[0.046]*bins #down/up
    pdf = [0.032]*bins

    sig_unc = {
            "lumi": lumi,
            "JES": JES,
            "pileup": pileup,
            "ggH_LHE_scale":higgsPt,
            "ggH_xsec":xsec,
            "ggH_pdf": pdf,
            "csc_DNN":dnn,
            "csc_time_spread":timespread,
    }
    if "dt" in category:
        sig_unc["dt_rpcBX"]=rpcBX
        sig_unc["dt_rpcMatch"]=rpcMatch
    return sig_unc
def make_datacard_2tag(outDataCardsDir,modelName,  signal_rate, normalization, bkg_rate, observation, bkg_unc, bkg_unc_name, sig_unc,signal_region, prefix):
    a,b,c,d = bkg_rate[0], bkg_rate[1], bkg_rate[2], bkg_rate[3]
    nSig = len(signal_rate.keys())
    text_file = open(outDataCardsDir+modelName+".txt", "w")
    text_file.write('# signal norm {0} \n'.format(normalization))

    text_file.write('imax {0} \n'.format(4))
    text_file.write('jmax {0} \n'.format(nSig))
    text_file.write('kmax * \n')
    text_file.write('shapes * * FAKE \n')


    text_file.write('--------------- \n')
    text_file.write('--------------- \n')
    text_file.write('bin \t chA \t chB \t chC \t chD \n')
    text_file.write('observation \t {0:6.2f} \t {1:6.2f} \t {2:6.2f} \t {3:6.2f} \n'.format(observation[0],observation[1],observation[2],observation[3]))
    text_file.write('------------------------------ \n')
    text_file.write('bin '+'\t chA ' * (1+nSig) + '\t chB ' * (1+nSig) +'\t chC '*(1+nSig) +'\t chD '*(1+nSig) +'\n')
    process_name = '\t '+ (' \t ').join(list(signal_rate.keys())) + '\t bkg '
    text_file.write('process ' + process_name * 4 + '\n')
    process_number = '\t '+ (' \t ').join(list((np.arange(nSig)*-1).astype(str))) + '\t 1'
    text_file.write('process ' + process_number * 4 + '\n')
    rate_string = 'rate'
    for i in range(4):# 4 bins
        for k,v in signal_rate.items():
            rate_string +='\t {0:e} '.format(v[i])
        rate_string += '\t 1 '
    text_file.write(rate_string+'\n')
    text_file.write('------------------------------ \n')

    text_file.write(prefix+'A   rateParam       chA     bkg      (@0*@2/@1)                    '+prefix+'B,'+prefix+'C,'+prefix+'D \n')
    if b == 0: text_file.write(prefix+'B   rateParam       chB     bkg     {0:.2f}        [0,{1:.2f}] \n'.format(b, c*7))
    else: text_file.write(prefix+'B   rateParam       chB     bkg     {0:.2f}        [0,{1:.2f}] \n'.format(b, b*7))
    text_file.write(prefix+'C   rateParam       chC     bkg     {0:.2f}        [0,{1:.2f}] \n'.format(c, c*7))
    if d == 0:text_file.write(prefix+'D   rateParam       chD     bkg     {0:.2f}        [0,{1:.2f}] \n'.format(d, c*7))
    else: text_file.write(prefix+'D   rateParam       chD     bkg     {0:.2f}        [0,{1:.2f}] \n'.format(d, d*7))


    for k,v in signal_rate.items():text_file.write('norm rateParam * {0} 1  \n'.format(k))
    sig_unc_name = list(sig_unc.keys())
    
    for k, v in sig_unc.items():
        
        if 'mc_stat' in k:
            print("here")
            for j, bin in enumerate(['A', 'B', 'C', 'D']):#bin
                before = (len(signal_rate.keys())+1)*j
                after = (len(signal_rate.keys())+1)*4-before-1
                if v[j] > 0.0: 
                    text_file.write(f'{k}_{bin} \t gmN ' +str(int(v[j]))+ '  '+'\t -  '*before + str(signal_rate['signal'][j]/int(v[j])) + '\t - '*after +'\n')

                        
        else:    
            unc_text = f'{k} \t lnN'
            if len(v)==4:#symmetric uncertainties
                for j in range(4):#bin
                    if v[j] == 0.0:unc_text += ' \t -'
                    else: unc_text += ' \t '+str(v[j]+1)
                    unc_text += '\t - '
            else:
                print(v)
                for j in range(4):#bin A, B, C, D
                    if  v[j] == 0.0 and v[j+4] == 0.0: unc_text += ' \t -'
                    else:unc_text += ' \t {0}/{1}'.format(1-v[j],1+v[j+4])
                    unc_text += '\t -'
            text_file.write(unc_text + ' \n')
            
    for i in range(len(bkg_unc_name)):
        bkg_unc_text = bkg_unc_name[i] + ' \t lnN ' + '\t - '*(4*nSig+3) + '\t ' + str(1+bkg_unc[i]) + ' \n'
        text_file.write(bkg_unc_text)
    

    text_file.close()
def readNorm(f_cscCard):
    f = open(f_cscCard,"r")
    norm = float(f.readline().split()[3])
    return norm


def HLT_CSC(eta,nstation,size):
    cond = (nstation == 1) & (np.abs(eta)<1.9) & (size >= 200)
    cond = cond | ((nstation == 1) & (np.abs(eta)>1.9) & (size >= 500))
    cond = cond | ((nstation > 1) & (np.abs(eta)<1.9) & (size > 100))
    cond = cond | ((nstation > 1) & (np.abs(eta)>1.9) & (size > 500))
    return cond

def L1_trg(cscClusterR, cscClusterZ, cscClusterSize):  
    first_in_ME11 = np.logical_and(np.logical_and(np.logical_and(cscClusterR>100, cscClusterR<275), np.abs(cscClusterZ)>580), np.abs(cscClusterZ)<632) 
    first_in_ME12 = np.logical_and(np.logical_and(np.logical_and(cscClusterR>275, cscClusterR<465), np.abs(cscClusterZ)>668), np.abs(cscClusterZ)<724)
    first_in_ME13 = np.logical_and(np.logical_and(np.logical_and(cscClusterR>505, cscClusterR<700), np.abs(cscClusterZ)>668), np.abs(cscClusterZ)<724)
    first_in_ME21 = np.logical_and(np.logical_and(np.logical_and(cscClusterR>139, cscClusterR<345), np.abs(cscClusterZ)>789), np.abs(cscClusterZ)<850)
    first_in_ME22 = np.logical_and(np.logical_and(np.logical_and(cscClusterR>357, cscClusterR<700), np.abs(cscClusterZ)>791), np.abs(cscClusterZ)<850)
    first_in_ME31 = np.logical_and(np.logical_and(np.logical_and(cscClusterR>160, cscClusterR<345), np.abs(cscClusterZ)>915), np.abs(cscClusterZ)<970)
    first_in_ME32 = np.logical_and(np.logical_and(np.logical_and(cscClusterR>357, cscClusterR<700), np.abs(cscClusterZ)>911), np.abs(cscClusterZ)<970)
    first_in_ME41 = np.logical_and(np.logical_and(np.logical_and(cscClusterR>178, cscClusterR<345), np.abs(cscClusterZ)>1002), np.abs(cscClusterZ)<1063)
    first_in_ME42 = np.logical_and(np.logical_and(np.logical_and(cscClusterR>357, cscClusterR<700), np.abs(cscClusterZ)>1002), np.abs(cscClusterZ)<1063)
    
    first_in_plateau_ME11 = np.logical_and(first_in_ME11, cscClusterSize>=500)
    first_in_plateau_ME21 = np.logical_and(first_in_ME21, cscClusterSize>=500)
    first_in_plateau_ME31 = np.logical_and(first_in_ME31, cscClusterSize>=500)
    first_in_plateau_ME41 = np.logical_and(first_in_ME41, cscClusterSize>=500)

    first_in_plateau_ME12 = np.logical_and(first_in_ME12, cscClusterSize>=200)
    first_in_plateau_ME13 = np.logical_and(first_in_ME13, cscClusterSize>=200)
    first_in_plateau_ME22 = np.logical_and(first_in_ME22, cscClusterSize>=200)
    first_in_plateau_ME32 = np.logical_and(first_in_ME32, cscClusterSize>=200)
    first_in_plateau_ME42 = np.logical_and(first_in_ME42, cscClusterSize>=200)
    
    first_in_plateau = first_in_plateau_ME11 | first_in_plateau_ME12 | first_in_plateau_ME13 | first_in_plateau_ME21 | first_in_plateau_ME22 | \
    first_in_plateau_ME31 | first_in_plateau_ME32 | first_in_plateau_ME41 | first_in_plateau_ME42
    return first_in_plateau

def deltaPhi( phi1,  phi2):
    dphi = phi1-phi2
    while np.count_nonzero(dphi > math.pi)>0:
        dphi[dphi > math.pi] -= 2*math.pi
    while np.count_nonzero(dphi< -math.pi)>0:
        dphi[dphi < -math.pi] += 2*math.pi
    return dphi
