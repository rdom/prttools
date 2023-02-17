import sys
import ROOT
import numpy as np

ROOT.gInterpreter.ProcessLine('#include "../../prttools/PrtTools.h"')
ROOT.gSystem.Load('../../prttools/PrtTools_cxx.so')

t = ROOT.PrtTools("beam_18215163732S.root")

stat = int(sys.argv[1])
ne_train = stat


x_train = np.zeros((ne_train,8,64))
y_train = np.zeros((ne_train,1))

x_test = np.zeros((4000,8,64))
y_test = np.zeros((4000,1))


while t.next() and t.i() < ne_train + 4000 :
    i = t.i()
    for hit in t.event().getHits() :
        # tof = e.getTof()
        # tofPi = fEvent->getTofPi()
        # tofP = fEvent->getTofP()
        time = hit.getLeadTime()
        if time > 30 :
            continue
        
        if i < ne_train:
            x_train[i,hit.getPmt(),hit.getPixel()-1] = 1 #time
            y_train[i] = t.pid()
        else :
            x_test[i-ne_train,hit.getPmt(),hit.getPixel()-1] = 1 #time
            y_test[i-ne_train] = t.pid()
        

nid = "data_stat_" + str(ne_train);
np.savez(nid, x_train= x_train, y_train=y_train, x_test=x_test, y_test= y_test)
