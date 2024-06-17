#!/cvmfs/cms.cern.ch/slc7_amd64_gcc700/cms/cmssw/CMSSW_10_2_13/external/slc7_amd64_gcc700/bin/python3
import os
import sys
import ROOT as rt
from random import shuffle
from time import time
import pickle
import numpy as np

datacards = [dc for dc in sys.argv[1:] if "datacard" in dc and "low" in dc and ".txt" in dc]
ctaus, mets, limits, limitErrs = [], [], [], []
for card in datacards:
    ct = card.split("_")[2]
    if ct not in ctaus:
        ctaus.append(ct)
    # ct = int(card.split("_")[1][2:].lstrip("0"))
    # if ct not in ctaus:
    #     ctaus.append(ct)
ctaus = sorted(ctaus, key=lambda x: abs(np.log10(int(x[2:].lstrip("0"))/999.99)))

channels = ("datacard", "csccsc", "cscdt", "2022", "2023", "csccsc_2022", "cscdt_2022", "csccsc_2023", "cscdt_2023")
limits = {
    channel: {
        k: {t: [] for t in ("ctau", "2.5", "16.0", "50.0", "84.0", "97.5", "limit")} for k in ("low", "high", "comb")
    }
    for channel in channels
}


def combine_cards(cards):
    #     # print('SKIPPING LOW/HIGH CARD COMBINATION')
    norms = []
    for c in cards:
        with open(c, "r") as fdc:
            norm = fdc.readline().rstrip().split(" ")
            norm = float(norm[-1])
            norms.append(norm)
    norm = 1/sum([1/norm for norm in norms if norm>0])
    
    if len(cards) > 1:
        card = "card_out.txt"
        print(cards)
        combine_cards_cmd = "combineCards.py "+" ".join(["Name"+str(i+1)+"="+dc for i, dc in enumerate(cards)])+" > "+card
        print(combine_cards_cmd)
        os.system(combine_cards_cmd)
    elif len(cards) == 1:
        card = cards[0]
    else:
        raise ValueError("cards empty")
    


    combine_cmd = "combine -M AsymptoticLimits --freezeParam norm --setParameters norm="+str(norm)+" "+card
    print(combine_cmd)
    os.system(combine_cmd)

    tfile = rt.TFile.Open("higgsCombineTest.AsymptoticLimits.mH120.root")
    limit_tr = tfile.limit

    # Multiply by norm to fix issues
    lms = [norm*lm for lm in [ ci.limit for ci in limit_tr]]
    lmEs = [norm*lmE for lmE in [ ci.limitErr for ci in limit_tr]]

    return [lms, lmEs]


time_str, tstart = "{:0>2.0f}:{:0>2.0f}:{:0>2.0f}", time()
for channel in channels:#limits.keys():
    channel_cards = [dc for dc in datacards if all(["_"+c in dc for c in channel.split("_")])]
    print(channel, len(channel_cards))
    for ct in ctaus:
        # LOW
        l_limits = combine_cards([cc for cc in channel_cards if ct in cc])

        # HIGH
        h_limits = combine_cards([cc.replace("_low", "_high") for cc in channel_cards if ct in cc])

        # COMB
        c_limits = combine_cards(
            [cc for cc in channel_cards if ct in cc]
            + [cc.replace("_low", "_high") for cc in channel_cards if ct in cc]
        )

        limits[channel]["low"]["ctau"].append(int(ct[2:].lstrip("0")))
        limits[channel]["high"]["ctau"].append(int(ct[2:].lstrip("0")))
        limits[channel]["comb"]["ctau"].append(int(ct[2:].lstrip("0")))
        for il, ll in enumerate(["2.5", "16.0", "50.0", "84.0", "97.5", "limit"]):
            limits[channel]["low"][ll].append(l_limits[0][il] if len(l_limits[0])>1 else 1)
            limits[channel]["high"][ll].append(h_limits[0][il] if len(h_limits[0])>1 else 1)
            limits[channel]["comb"][ll].append(c_limits[0][il] if len(c_limits[0])>1 else 1)

        with open('limits.pkl', 'wb') as fld:
            pickle.dump(limits, fld)

###################

# datacards = sys.argv[1:]
# ctaus, mets, limits, limitErrs = [], [], [], []

# limits = {
#     tag : {
#         k : {
#             t : [] for t in ('ctau','2.5','16.0','50.0','84.0','97.5','limit')
#         } for k in ('low','high','comb')
#     } for tag in ('csccsc','cscdt')
# }


# for dc in datacards:
#     if "datacard_" not in dc or "_low.txt" not in dc:
#         datacards.remove(dc)

# print(datacards)
# shuffle(datacards)
# # datacards = sorted(datacards, key=lambda x: abs(int(x.split('_')[2][2:])-1000))
# datacards = sorted(datacards, key=lambda x: abs(np.log10(int(x.split('_')[2][2:])/999.9)))

# time_str, tstart = "{:0>2.0f}:{:0>2.0f}:{:0>2.0f}", time()
# for idc, datacard in enumerate(datacards):
#     # if idc == 3:
#     #     break

#     # datacard_ct00100000_csccsc_low.txt
#     dc_pars = datacard.split("/")[-1][:-4].split("_")
#     ct = dc_pars[1][2:].lstrip("0")
#     tag = dc_pars[2]
#     region = dc_pars[3]
#     ct = ct if len(ct) else "0"

#     print('')
#     print(80*'=')
#     print(80*'=')
#     print(str(idc+1)+"/"+str(len(datacards))+" ("+str((idc+1)/len(datacards)*100)+"%)")
#     print(datacard)
#     tdone = (time() - tstart)
#     tleft = (len(datacards)-idc)/idc * tdone if idc else 0
#     print("  Time Elapsed (HH:MM:SS): " + time_str.format((tdone//60//60)%60, (tdone//60)%60, tdone%60))
#     print("  Time Left (HH:MM:SS):    " + time_str.format((tleft//60//60)%60, (tleft//60)%60, tleft%60))
#     print(40*'- ')
#     # *** #
#     datacard_low = datacard
#     datacard_high = datacard.replace("low","high")
#     datacard_comb = datacard.replace("low","comb")

#     # print('SKIPPING LOW/HIGH CARD COMBINATION')
#     combine_cards_cmd = "combineCards.py Name1="+datacard_low+" Name2="+datacard_high+" > "+datacard_comb
#     print(combine_cards_cmd)
#     os.system(combine_cards_cmd)

#     # *** #
#     norms = []
#     for card_id in ('low','high','comb'):
#         if card_id == 'all':
#             card = datacard_all
#         elif card_id == 'low':
#             card = datacard_low
#         elif card_id == 'high':
#             card = datacard_high
#         elif card_id == 'comb':
#             card = datacard_comb
#         else:
#             raise ValueError('oops')

#         # print('FINDING LIMITS FOR LOW MET CATEGORY ONLY')
#         # if 'low' in card_id:
#         if 'comb' not in card_id:
#             with open(card, "r") as fdc:
#                 norm = fdc.readline().rstrip().split(" ")
#                 norm = float(norm[-1])
#                 norms.append(norm)
#         else:
#             # norm = min(norms[-2:]) # sum(norms[-2:])/2
#             norm = 1/(1/norms[-1]+1/norms[-2])

#         combine_cmd = "combine -M AsymptoticLimits --freezeParam norm --setParameters norm="+str(norm)+" "+card
#         print(combine_cmd)
#         os.system(combine_cmd)

#         tfile = rt.TFile.Open("higgsCombineTest.AsymptoticLimits.mH120.root")
#         limit_tr = tfile.limit

#         # Multiply by norm to fix issues
#         lms = [norm*lm for lm in [ ci.limit for ci in limit_tr]]
#         lmEs = [norm*lmE for lmE in [ ci.limitErr for ci in limit_tr]]
#         limits[tag][card_id]['ctau'].append(int(ct))
#         limits[tag][card_id]['2.5'].append(lms[0])
#         limits[tag][card_id]['16.0'].append(lms[1])
#         limits[tag][card_id]['50.0'].append(lms[2])
#         limits[tag][card_id]['84.0'].append(lms[3])
#         limits[tag][card_id]['97.5'].append(lms[4])
#         limits[tag][card_id]['limit'].append(lms[5])

#         tfile.Close()
#         # else:
#         #     lms = [1 for lm in range(6)]
#         #     lmEs = [999 for lmE in range(6)]
#         #     limits[tag][card_id]['ctau'].append(int(ct))
#         #     limits[tag][card_id]['2.5'].append(lms[0])
#         #     limits[tag][card_id]['16.0'].append(lms[1])
#         #     limits[tag][card_id]['50.0'].append(lms[2])
#         #     limits[tag][card_id]['84.0'].append(lms[3])
#         #     limits[tag][card_id]['97.5'].append(lms[4])
#         #     limits[tag][card_id]['limit'].append(lms[5])

#     with open('limits.pkl', 'wb') as fld:
#         pickle.dump(limits, fld)
