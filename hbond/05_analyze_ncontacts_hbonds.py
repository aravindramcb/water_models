#!/usr/bin/env python3
# -*- coding:utf-8 -*-
#
#  Copyright 2021 Carlos Eduardo Sequeiros Borja <carseq@amu.edu.pl>
#  

import os
import pickle
import numpy as np
from sys import argv
import matplotlib.pyplot as plt
from collections import defaultdict
from dataclasses import dataclass, field

@dataclass(order=True)
class TEvent:
    tid: int = field(compare=True, hash=True)
    watid: int = 0
    frame: int = -1
    radius: float = 999.999
    dist2lig: float = 999.999
    fraction: float = -1.0
    event_type: str = ""
    supercluster: int = -1
    caver_cluster: int = -1

def normalized_plot_radii_by_SC(data, outfolder, tag="", frac_thr=0.7):
    mds = list(data.keys())
    mds.sort()
    
    sc_mds = defaultdict(int)
    for md in mds:
        _scs = set()
        for event in data[md]:
            if event.fraction >= frac_thr:
                _scs.add(event.supercluster)
        for s in _scs:
            sc_mds[s] += 1
    
    radii = defaultdict(list)
    for md in mds:
        for event in data[md]:
            if (event.fraction >= frac_thr) and (sc_mds[event.supercluster] >= 5):
                radii[event.supercluster].append(event.radius)
    
    scs = list(radii.keys())
    scs.sort()
    radii_hists = {}
    for sc in scs:
        radii_hists[sc] = np.histogram(radii[sc], range=(0.7, 3.7), bins=30)[0].astype(float)
        _, radii_xranges = np.histogram(radii[sc], range=(0.7, 3.7), bins=30)
    x_offset = (radii_xranges[1] - radii_xranges[0]) / 2.0
    hist_sums = np.zeros(30)
    for sc in scs:
        hist_sums += radii_hists[sc]
    for sc in scs:
        radii_hists[sc] = np.nan_to_num(radii_hists[sc]/hist_sums, nan=0.0)
    
    bottoms = np.zeros((len(scs), radii_hists[scs[0]].shape[0]))
    for i in range(1, len(scs)):
        for j in range(i, len(scs)):
            bottoms[j,:] += radii_hists[scs[i-1]]
    
    fig, ax = plt.subplots()
    for i in range(len(scs)):
        ax.bar(radii_xranges[:-1]+x_offset, height=radii_hists[scs[i]], width=0.08,
               bottom=bottoms[i], label="SC_{}".format(scs[i]))
    ax.set_xlim(radii_xranges[0]-0.01, radii_xranges[-1]+0.01)
    ax.set_ylim(0, 1)
    ax.set_xticks(radii_xranges)
    ax.set_xticklabels(labels=[str(np.round(x,1)) for x in radii_xranges],
                       fontdict={"fontsize":6})
    ax.legend(loc="upper right", fontsize=6)
    plt.ylabel("Fraction of events")
    plt.xlabel(r"Radii [$\AA$]")
    plt.title("Radii of closest sphere to water transport by SC in 5-MDs or more")
    plt.savefig(os.path.join(outfolder, tag+"radii_by_SC_normalized_md5.png"), dpi=300, format="png")

def normalized_plot_radii_by_HB(data,hbond_outfolder, outfolder, tag="", frac_thr=0.7):
    def _get_hbonds(hbondfile):
        hbonds = 0
        with open(hbondfile, "r") as hfile:
            for line in hfile:
                if line.startswith("#"):
                    continue
                hbonds = int(line.strip().split()[2])
        return hbonds
    
    mds = list(data.keys())
    mds.sort()
    for md in mds:
        for event in data[md]:
            hbondfile = "{}/{}/f{:0>5}_e{:0>4}.txt".format(hbond_outfolder, md, event.frame+1, event.tid)
            event.hbonds = _get_hbonds(hbondfile)
    
    sc_mds = defaultdict(int)
    for md in mds:
        _scs = set()
        for event in data[md]:
            if event.fraction >= frac_thr:
                _scs.add(event.supercluster)
        for s in _scs:
            sc_mds[s] += 1
    
    e_hb = defaultdict(list)
    for md in mds:
        for event in data[md]:
            if (event.fraction >= frac_thr) and (sc_mds[event.supercluster] >= 1):
                e_hb[event.hbonds].append(event.radius)
    
    ehb_hists = {}
    hbonds = list(e_hb.keys())
    hbonds.sort()
    for hb in hbonds:
        ehb_hists[hb] = np.histogram(e_hb[hb], range=(0.9, 2.7), bins=12)[0].astype(float)
        _, hb_xranges = np.histogram(e_hb[hb], range=(0.9, 2.7), bins=12)
    x_offset = (hb_xranges[1] - hb_xranges[0]) / 2.0
    hist_sums = np.zeros(12)
    for hb in hbonds:
        hist_sums += ehb_hists[hb]
    for hb in hbonds:
        ehb_hists[hb] = np.nan_to_num(ehb_hists[hb]/hist_sums, nan=0.0)
    
    bottoms = np.zeros((len(hbonds), ehb_hists[hbonds[0]].shape[0]))
    for i in range(1, len(hbonds)):
        for j in range(i, len(hbonds)):
            bottoms[j,:] += ehb_hists[hbonds[i-1]]
    
    fig, ax = plt.subplots()
    for i in range(len(hbonds)):
        ax.bar(hb_xranges[:-1]+x_offset, height=ehb_hists[hbonds[i]], width=0.09,
               bottom=bottoms[i], label="H-bonds_{}".format(hbonds[i]))
    ax.set_xlim(hb_xranges[0]-0.01, hb_xranges[-1]+0.01)
    ax.set_ylim(0, 1)
    ax.set_xticks(hb_xranges)
    ax.set_xticklabels(labels=[str(np.round(x,1)) for x in hb_xranges],
                       fontdict={"fontsize":6})
    ax.legend(loc="upper right", fontsize=6)
    plt.ylabel("Fraction of events")
    plt.xlabel(r"Radii [$\AA$]")
    plt.title("Radii of closest sphere to water transport by H-bonds in 1-MDs or more")
    plt.savefig(os.path.join(outfolder, tag+"radii_by_HB_normalized_md1.png"), dpi=300, format="png")

def normalized_plot_radii_by_restype(data, outfolder, threshold=3.0, tag="", frac_thr=0.7):
    dhaa_res = {1:"ILE", 2:"GLY", 3:"THR", 4:"GLY", 5:"PHE", 6:"PRO", 7:"PHE", 8:"ASP", 9:"PRO", 10:"HID", 11:"TYR", 12:"VAL", 13:"GLU", 14:"VAL", 15:"LEU", 16:"GLY", 17:"GLU", 18:"ARG", 19:"MET", 20:"HID", 21:"TYR", 22:"VAL", 23:"ASP", 24:"VAL", 25:"GLY", 26:"PRO", 27:"ARG", 28:"ASP", 29:"GLY", 30:"THR", 31:"PRO", 32:"VAL", 33:"LEU", 34:"PHE", 35:"LEU", 36:"HID", 37:"GLY", 38:"ASN", 39:"PRO", 40:"THR", 41:"SER", 42:"SER", 43:"TYR", 44:"LEU", 45:"TRP", 46:"ARG", 47:"ASN", 48:"ILE", 49:"ILE", 50:"PRO", 51:"HIE", 52:"VAL", 53:"ALA", 54:"PRO", 55:"SER", 56:"HIE", 57:"ARG", 58:"CYS", 59:"ILE", 60:"ALA", 61:"PRO", 62:"ASP", 63:"LEU", 64:"ILE", 65:"GLY", 66:"MET", 67:"GLY", 68:"LYS", 69:"SER", 70:"ASP", 71:"LYS", 72:"PRO", 73:"ASP", 74:"LEU", 75:"ASP", 76:"TYR", 77:"PHE", 78:"PHE", 79:"ASP", 80:"ASP", 81:"HIE", 82:"VAL", 83:"ARG", 84:"TYR", 85:"LEU", 86:"ASP", 87:"ALA", 88:"PHE", 89:"ILE", 90:"GLU", 91:"ALA", 92:"LEU", 93:"GLY", 94:"LEU", 95:"GLU", 96:"GLU", 97:"VAL", 98:"VAL", 99:"LEU", 100:"VAL", 101:"ILE", 102:"HIE", 103:"ASP", 104:"TRP", 105:"GLY", 106:"SER", 107:"ALA", 108:"LEU", 109:"GLY", 110:"PHE", 111:"HIE", 112:"TRP", 113:"ALA", 114:"LYS", 115:"ARG", 116:"ASN", 117:"PRO", 118:"GLU", 119:"ARG", 120:"VAL", 121:"LYS", 122:"GLY", 123:"ILE", 124:"ALA", 125:"CYS", 126:"MET", 127:"GLU", 128:"PHE", 129:"ILE", 130:"ARG", 131:"PRO", 132:"ILE", 133:"PRO", 134:"THR", 135:"TRP", 136:"ASP", 137:"GLU", 138:"TRP", 139:"PRO", 140:"GLU", 141:"PHE", 142:"ALA", 143:"ARG", 144:"GLU", 145:"THR", 146:"PHE", 147:"GLN", 148:"ALA", 149:"PHE", 150:"ARG", 151:"THR", 152:"ALA", 153:"ASP", 154:"VAL", 155:"GLY", 156:"ARG", 157:"GLU", 158:"LEU", 159:"ILE", 160:"ILE", 161:"ASP", 162:"GLN", 163:"ASN", 164:"ALA", 165:"PHE", 166:"ILE", 167:"GLU", 168:"GLY", 169:"ALA", 170:"LEU", 171:"PRO", 172:"LYS", 173:"CYS", 174:"VAL", 175:"VAL", 176:"ARG", 177:"PRO", 178:"LEU", 179:"THR", 180:"GLU", 181:"VAL", 182:"GLU", 183:"MET", 184:"ASP", 185:"HIE", 186:"TYR", 187:"ARG", 188:"GLU", 189:"PRO", 190:"PHE", 191:"LEU", 192:"LYS", 193:"PRO", 194:"VAL", 195:"ASP", 196:"ARG", 197:"GLU", 198:"PRO", 199:"LEU", 200:"TRP", 201:"ARG", 202:"PHE", 203:"PRO", 204:"ASN", 205:"GLU", 206:"LEU", 207:"PRO", 208:"ILE", 209:"ALA", 210:"GLY", 211:"GLU", 212:"PRO", 213:"ALA", 214:"ASN", 215:"ILE", 216:"VAL", 217:"ALA", 218:"LEU", 219:"VAL", 220:"GLU", 221:"ALA", 222:"TYR", 223:"MET", 224:"ASN", 225:"TRP", 226:"LEU", 227:"HID", 228:"GLN", 229:"SER", 230:"PRO", 231:"VAL", 232:"PRO", 233:"LYS", 234:"LEU", 235:"LEU", 236:"PHE", 237:"TRP", 238:"GLY", 239:"THR", 240:"PRO", 241:"GLY", 242:"VAL", 243:"LEU", 244:"ILE", 245:"PRO", 246:"PRO", 247:"ALA", 248:"GLU", 249:"ALA", 250:"ALA", 251:"ARG", 252:"LEU", 253:"ALA", 254:"GLU", 255:"SER", 256:"LEU", 257:"PRO", 258:"ASN", 259:"CYS", 260:"LYS", 261:"THR", 262:"VAL", 263:"ASP", 264:"ILE", 265:"GLY", 266:"PRO", 267:"GLY", 268:"LEU", 269:"HIP", 270:"TYR", 271:"LEU", 272:"GLN", 273:"GLU", 274:"ASP", 275:"ASN", 276:"PRO", 277:"ASP", 278:"LEU", 279:"ILE", 280:"GLY", 281:"SER", 282:"GLU", 283:"ILE", 284:"ALA", 285:"ARG", 286:"TRP", 287:"LEU", 288:"PRO", 289:"ALA", 290:"LEU", 291:"HIE", 292:"HIE", 293:"HIP"}
    epx_res = {1:"LYS", 2:"ILE", 3:"GLU", 4:"HIE", 5:"LYS", 6:"MET", 7:"VAL", 8:"ALA", 9:"VAL", 10:"ASN", 11:"GLY", 12:"LEU", 13:"ASN", 14:"MET", 15:"HIP", 16:"LEU", 17:"ALA", 18:"GLU", 19:"LEU", 20:"GLY", 21:"GLU", 22:"GLY", 23:"PRO", 24:"THR", 25:"ILE", 26:"LEU", 27:"PHE", 28:"ILE", 29:"HID", 30:"GLY", 31:"PHE", 32:"PRO", 33:"GLU", 34:"LEU", 35:"TRP", 36:"TYR", 37:"SER", 38:"TRP", 39:"ARG", 40:"HIE", 41:"GLN", 42:"MET", 43:"VAL", 44:"TYR", 45:"LEU", 46:"ALA", 47:"GLU", 48:"ARG", 49:"GLY", 50:"TYR", 51:"ARG", 52:"ALA", 53:"VAL", 54:"ALA", 55:"PRO", 56:"ASP", 57:"LEU", 58:"ARG", 59:"GLY", 60:"TYR", 61:"GLY", 62:"ASP", 63:"THR", 64:"THR", 65:"GLY", 66:"ALA", 67:"PRO", 68:"LEU", 69:"ASN", 70:"ASP", 71:"PRO", 72:"SER", 73:"LYS", 74:"PHE", 75:"SER", 76:"ILE", 77:"LEU", 78:"HIE", 79:"LEU", 80:"VAL", 81:"GLY", 82:"ASP", 83:"VAL", 84:"VAL", 85:"ALA", 86:"LEU", 87:"LEU", 88:"GLU", 89:"ALA", 90:"ILE", 91:"ALA", 92:"PRO", 93:"ASN", 94:"GLU", 95:"GLU", 96:"LYS", 97:"VAL", 98:"PHE", 99:"VAL", 100:"VAL", 101:"ALA", 102:"HIE", 103:"ASP", 104:"TRP", 105:"GLY", 106:"ALA", 107:"LEU", 108:"ILE", 109:"ALA", 110:"TRP", 111:"HID", 112:"LEU", 113:"CYS", 114:"LEU", 115:"PHE", 116:"ARG", 117:"PRO", 118:"ASP", 119:"LYS", 120:"VAL", 121:"LYS", 122:"ALA", 123:"LEU", 124:"VAL", 125:"ASN", 126:"LEU", 127:"SER", 128:"VAL", 129:"HIE", 130:"PHE", 131:"SER", 132:"LYS", 133:"ARG", 134:"ASN", 135:"PRO", 136:"LYS", 137:"MET", 138:"ASN", 139:"VAL", 140:"VAL", 141:"GLU", 142:"GLY", 143:"LEU", 144:"LYS", 145:"ALA", 146:"ILE", 147:"TYR", 148:"GLY", 149:"GLU", 150:"ASP", 151:"HIE", 152:"TYR", 153:"ILE", 154:"SER", 155:"ARG", 156:"PHE", 157:"GLN", 158:"VAL", 159:"PRO", 160:"GLY", 161:"GLU", 162:"ILE", 163:"GLU", 164:"ALA", 165:"GLU", 166:"PHE", 167:"ALA", 168:"PRO", 169:"ILE", 170:"GLY", 171:"ALA", 172:"LYS", 173:"SER", 174:"VAL", 175:"LEU", 176:"LYS", 177:"LYS", 178:"ILE", 179:"LEU", 180:"THR", 181:"TYR", 182:"ARG", 183:"ASP", 184:"PRO", 185:"ALA", 186:"PRO", 187:"PHE", 188:"TYR", 189:"PHE", 190:"PRO", 191:"LYS", 192:"GLY", 193:"LYS", 194:"GLY", 195:"LEU", 196:"GLU", 197:"ALA", 198:"ILE", 199:"PRO", 200:"ASP", 201:"ALA", 202:"PRO", 203:"VAL", 204:"ALA", 205:"LEU", 206:"SER", 207:"SER", 208:"TRP", 209:"LEU", 210:"SER", 211:"GLU", 212:"GLU", 213:"GLU", 214:"LEU", 215:"ASP", 216:"TYR", 217:"TYR", 218:"ALA", 219:"ASN", 220:"LYS", 221:"PHE", 222:"GLU", 223:"GLN", 224:"THR", 225:"GLY", 226:"PHE", 227:"THR", 228:"GLY", 229:"ALA", 230:"VAL", 231:"ASN", 232:"TYR", 233:"TYR", 234:"ARG", 235:"ALA", 236:"LEU", 237:"PRO", 238:"ILE", 239:"ASN", 240:"TRP", 241:"GLU", 242:"LEU", 243:"THR", 244:"ALA", 245:"PRO", 246:"TRP", 247:"THR", 248:"GLY", 249:"ALA", 250:"GLN", 251:"VAL", 252:"LYS", 253:"VAL", 254:"PRO", 255:"THR", 256:"LYS", 257:"PHE", 258:"ILE", 259:"VAL", 260:"GLY", 261:"GLU", 262:"PHE", 263:"ASP", 264:"LEU", 265:"VAL", 266:"TYR", 267:"HIP", 268:"ILE", 269:"PRO", 270:"GLY", 271:"ALA", 272:"LYS", 273:"GLU", 274:"TYR", 275:"ILE", 276:"HIE", 277:"ASN", 278:"GLY", 279:"GLY", 280:"PHE", 281:"LYS", 282:"LYS", 283:"ASP", 284:"VAL", 285:"PRO", 286:"LEU", 287:"LEU", 288:"GLU", 289:"GLU", 290:"VAL", 291:"VAL", 292:"VAL", 293:"LEU", 294:"GLU", 295:"GLY", 296:"ALA", 297:"ALA", 298:"HIP", 299:"PHE", 300:"VAL", 301:"SER", 302:"GLN", 303:"GLU", 304:"ARG", 305:"PRO", 306:"HIE", 307:"GLU", 308:"ILE", 309:"SER", 310:"LYS", 311:"HIE", 312:"ILE", 313:"TYR", 314:"ASP", 315:"PHE", 316:"ILE", 317:"GLN", 318:"LYS", 319:"PHE"}
    lipase_res = {1:"ALA", 2:"PRO", 3:"THR", 4:"ALA", 5:"THR", 6:"LEU", 7:"ALA", 8:"ASN", 9:"GLY", 10:"ASP", 11:"THR", 12:"ILE", 13:"THR", 14:"GLY", 15:"LEU", 16:"ASN", 17:"ALA", 18:"ILE", 19:"ILE", 20:"ASN", 21:"GLU", 22:"ALA", 23:"PHE", 24:"LEU", 25:"GLY", 26:"ILE", 27:"PRO", 28:"PHE", 29:"ALA", 30:"GLU", 31:"PRO", 32:"PRO", 33:"VAL", 34:"GLY", 35:"ASN", 36:"LEU", 37:"ARG", 38:"PHE", 39:"LYS", 40:"ASP", 41:"PRO", 42:"VAL", 43:"PRO", 44:"TYR", 45:"SER", 46:"GLY", 47:"SER", 48:"LEU", 49:"ASP", 50:"GLY", 51:"GLN", 52:"LYS", 53:"PHE", 54:"THR", 55:"SER", 56:"TYR", 57:"GLY", 58:"PRO", 59:"SER", 60:"CYX", 61:"MET", 62:"GLN", 63:"GLN", 64:"ASN", 65:"PRO", 66:"GLU", 67:"GLY", 68:"THR", 69:"TYR", 70:"GLU", 71:"GLU", 72:"ASN", 73:"LEU", 74:"PRO", 75:"LYS", 76:"ALA", 77:"ALA", 78:"LEU", 79:"ASP", 80:"LEU", 81:"VAL", 82:"MET", 83:"GLN", 84:"SER", 85:"LYS", 86:"VAL", 87:"PHE", 88:"GLU", 89:"ALA", 90:"VAL", 91:"SER", 92:"PRO", 93:"SER", 94:"SER", 95:"GLU", 96:"ASP", 97:"CYX", 98:"LEU", 99:"THR", 100:"ILE", 101:"ASN", 102:"VAL", 103:"VAL", 104:"ARG", 105:"PRO", 106:"PRO", 107:"GLY", 108:"THR", 109:"LYS", 110:"ALA", 111:"GLY", 112:"ALA", 113:"ASN", 114:"LEU", 115:"PRO", 116:"VAL", 117:"MET", 118:"LEU", 119:"TRP", 120:"ILE", 121:"PHE", 122:"GLY", 123:"GLY", 124:"GLY", 125:"PHE", 126:"GLH", 127:"VAL", 128:"GLY", 129:"GLY", 130:"THR", 131:"SER", 132:"THR", 133:"PHE", 134:"PRO", 135:"PRO", 136:"ALA", 137:"GLN", 138:"MET", 139:"ILE", 140:"THR", 141:"LYS", 142:"SER", 143:"ILE", 144:"ALA", 145:"MET", 146:"GLY", 147:"LYS", 148:"PRO", 149:"ILE", 150:"ILE", 151:"HID", 152:"VAL", 153:"SER", 154:"VAL", 155:"ASN", 156:"TYR", 157:"ARG", 158:"VAL", 159:"SER", 160:"SER", 161:"TRP", 162:"GLY", 163:"PHE", 164:"LEU", 165:"ALA", 166:"GLY", 167:"ASP", 168:"GLU", 169:"ILE", 170:"LYS", 171:"ALA", 172:"GLU", 173:"GLY", 174:"SER", 175:"ALA", 176:"ASN", 177:"ALA", 178:"GLY", 179:"LEU", 180:"LYS", 181:"ASP", 182:"GLN", 183:"ARG", 184:"LEU", 185:"GLY", 186:"MET", 187:"GLN", 188:"TRP", 189:"VAL", 190:"ALA", 191:"ASP", 192:"ASN", 193:"ILE", 194:"ALA", 195:"ALA", 196:"PHE", 197:"GLY", 198:"GLY", 199:"ASP", 200:"PRO", 201:"THR", 202:"LYS", 203:"VAL", 204:"THR", 205:"ILE", 206:"PHE", 207:"GLY", 208:"GLU", 209:"SER", 210:"ALA", 211:"GLY", 212:"SER", 213:"MET", 214:"SER", 215:"VAL", 216:"MET", 217:"CYS", 218:"HIE", 219:"ILE", 220:"LEU", 221:"TRP", 222:"ASN", 223:"ASP", 224:"GLY", 225:"ASP", 226:"ASN", 227:"THR", 228:"TYR", 229:"LYS", 230:"GLY", 231:"LYS", 232:"PRO", 233:"LEU", 234:"PHE", 235:"ARG", 236:"ALA", 237:"GLY", 238:"ILE", 239:"MET", 240:"GLN", 241:"SER", 242:"GLY", 243:"ALA", 244:"MET", 245:"VAL", 246:"PRO", 247:"SER", 248:"ASP", 249:"ALA", 250:"VAL", 251:"ASP", 252:"GLY", 253:"ILE", 254:"TYR", 255:"GLY", 256:"ASN", 257:"GLU", 258:"ILE", 259:"PHE", 260:"ASP", 261:"LEU", 262:"LEU", 263:"ALA", 264:"SER", 265:"ASN", 266:"ALA", 267:"GLY", 268:"CYX", 269:"GLY", 270:"SER", 271:"ALA", 272:"SER", 273:"ASP", 274:"LYS", 275:"LEU", 276:"ALA", 277:"CYX", 278:"LEU", 279:"ARG", 280:"GLY", 281:"VAL", 282:"SER", 283:"SER", 284:"ASP", 285:"THR", 286:"LEU", 287:"GLH", 288:"ASP", 289:"ALA", 290:"THR", 291:"ASN", 292:"ASN", 293:"THR", 294:"PRO", 295:"GLY", 296:"PHE", 297:"LEU", 298:"ALA", 299:"TYR", 300:"SER", 301:"SER", 302:"LEU", 303:"ARG", 304:"LEU", 305:"SER", 306:"TYR", 307:"LEU", 308:"PRO", 309:"ARG", 310:"PRO", 311:"ASP", 312:"GLY", 313:"VAL", 314:"ASN", 315:"ILE", 316:"THR", 317:"ASP", 318:"ASP", 319:"MET", 320:"TYR", 321:"ALA", 322:"LEU", 323:"VAL", 324:"ARG", 325:"GLU", 326:"GLY", 327:"LYS", 328:"TYR", 329:"ALA", 330:"ASN", 331:"ILE", 332:"PRO", 333:"VAL", 334:"ILE", 335:"ILE", 336:"GLY", 337:"ASP", 338:"GLN", 339:"ASN", 340:"ASP", 341:"GLU", 342:"GLY", 343:"THR", 344:"PHE", 345:"PHE", 346:"GLY", 347:"THR", 348:"SER", 349:"SER", 350:"LEU", 351:"ASN", 352:"VAL", 353:"THR", 354:"THR", 355:"ASP", 356:"ALA", 357:"GLN", 358:"ALA", 359:"ARG", 360:"GLU", 361:"TYR", 362:"PHE", 363:"LYS", 364:"GLN", 365:"SER", 366:"PHE", 367:"VAL", 368:"HIP", 369:"ALA", 370:"SER", 371:"ASP", 372:"ALA", 373:"GLU", 374:"ILE", 375:"ASP", 376:"THR", 377:"LEU", 378:"MET", 379:"THR", 380:"ALA", 381:"TYR", 382:"PRO", 383:"GLY", 384:"ASP", 385:"ILE", 386:"THR", 387:"GLN", 388:"GLY", 389:"SER", 390:"PRO", 391:"PHE", 392:"ASP", 393:"THR", 394:"GLY", 395:"ILE", 396:"LEU", 397:"ASN", 398:"ALA", 399:"LEU", 400:"THR", 401:"PRO", 402:"GLN", 403:"PHE", 404:"LYS", 405:"ARG", 406:"ILE", 407:"SER", 408:"ALA", 409:"VAL", 410:"LEU", 411:"GLY", 412:"ASP", 413:"LEU", 414:"GLY", 415:"PHE", 416:"THR", 417:"LEU", 418:"ALA", 419:"ARG", 420:"ARG", 421:"TYR", 422:"PHE", 423:"LEU", 424:"ASN", 425:"HID", 426:"TYR", 427:"THR", 428:"GLY", 429:"GLY", 430:"THR", 431:"LYS", 432:"TYR", 433:"SER", 434:"PHE", 435:"LEU", 436:"SER", 437:"LYS", 438:"GLN", 439:"LEU", 440:"SER", 441:"GLY", 442:"LEU", 443:"PRO", 444:"VAL", 445:"LEU", 446:"GLY", 447:"THR", 448:"PHE", 449:"HIP", 450:"SER", 451:"ASN", 452:"ASP", 453:"ILE", 454:"VAL", 455:"PHE", 456:"GLN", 457:"ASH", 458:"TYR", 459:"LEU", 460:"LEU", 461:"GLY", 462:"SER", 463:"GLY", 464:"SER", 465:"LEU", 466:"ILE", 467:"TYR", 468:"ASN", 469:"ASN", 470:"ALA", 471:"PHE", 472:"ILE", 473:"ALA", 474:"PHE", 475:"ALA", 476:"THR", 477:"ASP", 478:"LEU", 479:"ASP", 480:"PRO", 481:"ASN", 482:"THR", 483:"ALA", 484:"GLY", 485:"LEU", 486:"LEU", 487:"VAL", 488:"LYS", 489:"TRP", 490:"PRO", 491:"GLU", 492:"TYR", 493:"THR", 494:"SER", 495:"SER", 496:"SER", 497:"GLN", 498:"SER", 499:"GLY", 500:"ASN", 501:"ASN", 502:"LEU", 503:"MET", 504:"MET", 505:"ILE", 506:"ASN", 507:"ALA", 508:"LEU", 509:"GLY", 510:"LEU", 511:"TYR", 512:"THR", 513:"GLY", 514:"LYS", 515:"ASP", 516:"ASN", 517:"PHE", 518:"ARG", 519:"THR", 520:"ALA", 521:"GLY", 522:"TYR", 523:"ASP", 524:"ALA", 525:"LEU", 526:"PHE", 527:"SER", 528:"ASN", 529:"PRO", 530:"PRO", 531:"SER", 532:"PHE", 533:"PHE", 534:"VAL"}
    hepx_wt_res = {1:"THR", 2:"SER", 3:"CYS", 4:"ASN", 5:"PRO", 6:"SER", 7:"ASP", 8:"MET", 9:"SER", 10:"HIE", 11:"GLY", 12:"TYR", 13:"VAL", 14:"THR", 15:"VAL", 16:"LYS", 17:"PRO", 18:"ARG", 19:"VAL", 20:"ARG", 21:"LEU", 22:"HIP", 23:"PHE", 24:"VAL", 25:"GLU", 26:"LEU", 27:"GLY", 28:"SER", 29:"GLY", 30:"PRO", 31:"ALA", 32:"VAL", 33:"CYS", 34:"LEU", 35:"CYS", 36:"HID", 37:"GLY", 38:"PHE", 39:"PRO", 40:"GLU", 41:"SER", 42:"TRP", 43:"TYR", 44:"SER", 45:"TRP", 46:"ARG", 47:"TYR", 48:"GLN", 49:"ILE", 50:"PRO", 51:"ALA", 52:"LEU", 53:"ALA", 54:"GLN", 55:"ALA", 56:"GLY", 57:"TYR", 58:"ARG", 59:"VAL", 60:"LEU", 61:"ALA", 62:"MET", 63:"ASP", 64:"MET", 65:"LYS", 66:"GLY", 67:"TYR", 68:"GLY", 69:"GLU", 70:"SER", 71:"SER", 72:"ALA", 73:"PRO", 74:"PRO", 75:"GLU", 76:"ILE", 77:"GLU", 78:"GLU", 79:"TYR", 80:"CYS", 81:"MET", 82:"GLU", 83:"VAL", 84:"LEU", 85:"CYS", 86:"LYS", 87:"GLU", 88:"MET", 89:"VAL", 90:"THR", 91:"PHE", 92:"LEU", 93:"ASP", 94:"LYS", 95:"LEU", 96:"GLY", 97:"LEU", 98:"SER", 99:"GLN", 100:"ALA", 101:"VAL", 102:"PHE", 103:"ILE", 104:"GLY", 105:"HIE", 106:"ASP", 107:"TRP", 108:"GLY", 109:"GLY", 110:"MET", 111:"LEU", 112:"VAL", 113:"TRP", 114:"TYR", 115:"MET", 116:"ALA", 117:"LEU", 118:"PHE", 119:"TYR", 120:"PRO", 121:"GLU", 122:"ARG", 123:"VAL", 124:"ARG", 125:"ALA", 126:"VAL", 127:"ALA", 128:"SER", 129:"LEU", 130:"ASN", 131:"THR", 132:"PRO", 133:"PHE", 134:"ILE", 135:"PRO", 136:"ALA", 137:"ASN", 138:"PRO", 139:"ASN", 140:"MET", 141:"SER", 142:"PRO", 143:"LEU", 144:"GLU", 145:"SER", 146:"ILE", 147:"LYS", 148:"ALA", 149:"ASN", 150:"PRO", 151:"VAL", 152:"PHE", 153:"ASP", 154:"TYR", 155:"GLN", 156:"LEU", 157:"TYR", 158:"PHE", 159:"GLN", 160:"GLU", 161:"PRO", 162:"GLY", 163:"VAL", 164:"ALA", 165:"GLU", 166:"ALA", 167:"GLU", 168:"LEU", 169:"GLU", 170:"GLN", 171:"ASN", 172:"LEU", 173:"SER", 174:"ARG", 175:"THR", 176:"PHE", 177:"LYS", 178:"SER", 179:"LEU", 180:"PHE", 181:"ARG", 182:"ALA", 183:"SER", 184:"ASP", 185:"GLU", 186:"SER", 187:"VAL", 188:"LEU", 189:"SER", 190:"MET", 191:"HIE", 192:"LYS", 193:"VAL", 194:"CYS", 195:"GLU", 196:"ALA", 197:"GLY", 198:"GLY", 199:"LEU", 200:"PHE", 201:"VAL", 202:"ASN", 203:"SER", 204:"PRO", 205:"GLU", 206:"GLU", 207:"PRO", 208:"SER", 209:"LEU", 210:"SER", 211:"ARG", 212:"MET", 213:"VAL", 214:"THR", 215:"GLU", 216:"GLU", 217:"GLU", 218:"ILE", 219:"GLN", 220:"PHE", 221:"TYR", 222:"VAL", 223:"GLN", 224:"GLN", 225:"PHE", 226:"LYS", 227:"LYS", 228:"SER", 229:"GLY", 230:"PHE", 231:"ARG", 232:"GLY", 233:"PRO", 234:"LEU", 235:"ASN", 236:"TRP", 237:"TYR", 238:"ARG", 239:"ASN", 240:"MET", 241:"GLU", 242:"ARG", 243:"ASN", 244:"TRP", 245:"LYS", 246:"TRP", 247:"ALA", 248:"CYS", 249:"LYS", 250:"SER", 251:"LEU", 252:"GLY", 253:"ARG", 254:"LYS", 255:"ILE", 256:"LEU", 257:"ILE", 258:"PRO", 259:"ALA", 260:"LEU", 261:"MET", 262:"VAL", 263:"THR", 264:"ALA", 265:"GLU", 266:"LYS", 267:"ASP", 268:"PHE", 269:"VAL", 270:"LEU", 271:"VAL", 272:"PRO", 273:"GLN", 274:"MET", 275:"SER", 276:"GLN", 277:"HIE", 278:"MET", 279:"GLU", 280:"ASP", 281:"TRP", 282:"ILE", 283:"PRO", 284:"HIE", 285:"LEU", 286:"LYS", 287:"ARG", 288:"GLY", 289:"HIE", 290:"ILE", 291:"GLU", 292:"ASP", 293:"CYS", 294:"GLY", 295:"HIP", 296:"TRP", 297:"THR", 298:"GLN", 299:"MET", 300:"ASP", 301:"LYS", 302:"PRO", 303:"THR", 304:"GLU", 305:"VAL", 306:"ASN", 307:"GLN", 308:"ILE", 309:"LEU", 310:"ILE", 311:"LYS", 312:"TRP", 313:"LEU", 314:"ASP", 315:"SER", 316:"ASP"}
    hepx_mut_res = {1:"THR", 2:"SER", 3:"CYS", 4:"ASN", 5:"PRO", 6:"SER", 7:"ASP", 8:"MET", 9:"SER", 10:"HIE", 11:"GLY", 12:"TYR", 13:"VAL", 14:"THR", 15:"VAL", 16:"LYS", 17:"PRO", 18:"ARG", 19:"VAL", 20:"ARG", 21:"LEU", 22:"HIP", 23:"PHE", 24:"VAL", 25:"GLU", 26:"LEU", 27:"GLY", 28:"SER", 29:"GLY", 30:"PRO", 31:"ALA", 32:"VAL", 33:"CYS", 34:"LEU", 35:"CYS", 36:"HID", 37:"GLY", 38:"PHE", 39:"PRO", 40:"GLU", 41:"SER", 42:"TRP", 43:"TYR", 44:"SER", 45:"TRP", 46:"ARG", 47:"TYR", 48:"GLN", 49:"ILE", 50:"PRO", 51:"ALA", 52:"LEU", 53:"ALA", 54:"GLN", 55:"ALA", 56:"GLY", 57:"TYR", 58:"ARG", 59:"VAL", 60:"LEU", 61:"ALA", 62:"MET", 63:"ASP", 64:"MET", 65:"LYS", 66:"GLY", 67:"TYR", 68:"GLY", 69:"GLU", 70:"SER", 71:"SER", 72:"ALA", 73:"PRO", 74:"PRO", 75:"GLU", 76:"ILE", 77:"GLU", 78:"GLU", 79:"TYR", 80:"CYS", 81:"MET", 82:"GLU", 83:"VAL", 84:"LEU", 85:"CYS", 86:"LYS", 87:"GLU", 88:"MET", 89:"VAL", 90:"THR", 91:"PHE", 92:"LEU", 93:"ASP", 94:"LYS", 95:"LEU", 96:"GLY", 97:"LEU", 98:"SER", 99:"GLN", 100:"ALA", 101:"VAL", 102:"PHE", 103:"ILE", 104:"GLY", 105:"HIE", 106:"ASP", 107:"TRP", 108:"GLY", 109:"GLY", 110:"MET", 111:"LEU", 112:"VAL", 113:"TRP", 114:"TYR", 115:"MET", 116:"ALA", 117:"LEU", 118:"PHE", 119:"TYR", 120:"PRO", 121:"GLU", 122:"ARG", 123:"VAL", 124:"ARG", 125:"ALA", 126:"VAL", 127:"ALA", 128:"SER", 129:"LEU", 130:"ASN", 131:"THR", 132:"PRO", 133:"PHE", 134:"ILE", 135:"PRO", 136:"ALA", 137:"ASN", 138:"PRO", 139:"ASN", 140:"MET", 141:"SER", 142:"PRO", 143:"LEU", 144:"GLU", 145:"SER", 146:"ILE", 147:"LYS", 148:"ALA", 149:"ASN", 150:"PRO", 151:"VAL", 152:"PHE", 153:"ASP", 154:"TYR", 155:"GLN", 156:"LEU", 157:"TYR", 158:"PHE", 159:"GLN", 160:"GLU", 161:"PRO", 162:"GLY", 163:"VAL", 164:"ALA", 165:"GLU", 166:"ALA", 167:"GLU", 168:"LEU", 169:"GLU", 170:"GLN", 171:"ASN", 172:"LEU", 173:"SER", 174:"ARG", 175:"THR", 176:"PHE", 177:"LYS", 178:"SER", 179:"LEU", 180:"PHE", 181:"ARG", 182:"ALA", 183:"SER", 184:"ASP", 185:"GLU", 186:"SER", 187:"VAL", 188:"LEU", 189:"SER", 190:"MET", 191:"HIE", 192:"LYS", 193:"VAL", 194:"CYS", 195:"GLU", 196:"ALA", 197:"GLY", 198:"GLY", 199:"LEU", 200:"PHE", 201:"VAL", 202:"ASN", 203:"SER", 204:"PRO", 205:"GLU", 206:"GLU", 207:"PRO", 208:"SER", 209:"LEU", 210:"SER", 211:"ARG", 212:"MET", 213:"VAL", 214:"THR", 215:"GLU", 216:"GLU", 217:"GLU", 218:"ILE", 219:"GLN", 220:"PHE", 221:"TYR", 222:"VAL", 223:"GLN", 224:"GLN", 225:"PHE", 226:"LYS", 227:"LYS", 228:"SER", 229:"GLY", 230:"PHE", 231:"ARG", 232:"GLY", 233:"PRO", 234:"LEU", 235:"ASN", 236:"TRP", 237:"TYR", 238:"ARG", 239:"ASN", 240:"MET", 241:"GLY", 242:"ARG", 243:"ASN", 244:"TRP", 245:"LYS", 246:"TRP", 247:"ALA", 248:"CYS", 249:"LYS", 250:"SER", 251:"LEU", 252:"GLY", 253:"ARG", 254:"LYS", 255:"ILE", 256:"LEU", 257:"ILE", 258:"PRO", 259:"ALA", 260:"LEU", 261:"MET", 262:"VAL", 263:"THR", 264:"ALA", 265:"GLU", 266:"LYS", 267:"ASP", 268:"PHE", 269:"VAL", 270:"LEU", 271:"VAL", 272:"PRO", 273:"GLN", 274:"MET", 275:"SER", 276:"GLN", 277:"HIE", 278:"MET", 279:"GLU", 280:"ASP", 281:"TRP", 282:"ILE", 283:"PRO", 284:"HIE", 285:"LEU", 286:"LYS", 287:"ARG", 288:"GLY", 289:"HIE", 290:"ILE", 291:"GLU", 292:"ASP", 293:"CYS", 294:"GLY", 295:"HIP", 296:"TRP", 297:"THR", 298:"GLN", 299:"MET", 300:"ASP", 301:"LYS", 302:"PRO", 303:"THR", 304:"GLU", 305:"VAL", 306:"ASN", 307:"GLN", 308:"ILE", 309:"LEU", 310:"ILE", 311:"LYS", 312:"TRP", 313:"LEU", 314:"ASP", 315:"SER", 316:"ASP"}
    aa_neg = {"GLU", "GLH", "ASP", "ASH"}
    aa_pos = {"ARG", "HID", "HIE", "HIP", "LYS"}
    aa_polar = {"ASN", "CYS", "CYX", "GLN", "SER", "THR", "TYR"}
    aa_nonp = {"ALA", "GLY", "ILE", "LEU", "MET", "PHE", "PRO", "TRP", "VAL"}
    bbone_atoms = {"N", "CA", "C", "O", "NH", "HA", "H"}
    
    mds = list(data.keys())
    mds.sort()
    
    sc_mds = defaultdict(int)
    for md in mds:
        _scs = set()
        for event in data[md]:
            if event.fraction >= frac_thr:
                _scs.add(event.supercluster)
        for s in _scs:
            sc_mds[s] += 1
    
    residues = defaultdict(int)
    radii_by_type = defaultdict(list)
    prot = lipase_res
    for md in mds:
        for event in data[md]:
            _resids = set()
            contactfile = "{}/ncontacts/{}/f{:0>5}_e{:0>4}.txt".format(outfolder, md, event.frame+1, event.tid)
            if sc_mds[event.supercluster] >= 5:
                with open(contactfile, "r") as fin:
                    [fin.readline() for i in range(3)]
                    for line in fin:
                        chunks = line.strip().split()
                        dist = float(chunks[4])
                        _res = chunks[1][1:].split(":")[-1]
                        _res, _atom = _res.split("@")
                        if dist <= threshold:
                            r = int(_res)
                            if r not in _resids:
                                _resids.add(r)
                                if _atom in bbone_atoms:
                                    radii_by_type["bbone"].append(event.radius)
                                elif prot[r] in aa_neg:
                                    radii_by_type["aa_neg"].append(event.radius)
                                elif prot[r] in aa_pos:
                                    radii_by_type["aa_pos"].append(event.radius)
                                elif prot[r] in aa_polar:
                                    radii_by_type["aa_polar"].append(event.radius)
                                elif prot[r] in aa_nonp:
                                    radii_by_type["aa_nonp"].append(event.radius)
                                else:
                                    print("Residue not in groups:", prot[r])
    
    ctc_hists = {}
    ctc_types = ["bbone", "aa_nonp", "aa_polar", "aa_pos", "aa_neg"]
    colors = ["black", "tab:gray", "tab:olive", "tab:blue", "tab:red"]
    for ctc_type in ctc_types:
        ctc_hists[ctc_type] = np.histogram(radii_by_type[ctc_type], range=(0.7, 3.7), bins=30)[0].astype(float)
        _, ctc_xranges = np.histogram(radii_by_type[ctc_type], range=(0.7, 3.7), bins=30)
    x_offset = (ctc_xranges[1] - ctc_xranges[0]) / 2.0
    
    hist_sums = np.zeros(30)
    for ctc_type in ctc_types:
        hist_sums += ctc_hists[ctc_type]
    for ctc_type in ctc_types:
        ctc_hists[ctc_type] = np.nan_to_num(ctc_hists[ctc_type]/hist_sums, nan=0.0)
    
    bottoms = np.zeros((len(ctc_types), ctc_hists[ctc_types[0]].shape[0]))
    for i in range(1, len(ctc_types)):
        for j in range(i, len(ctc_types)):
            bottoms[j,:] += ctc_hists[ctc_types[i-1]]
    
    fg, ax = plt.subplots()
    for i in range(len(ctc_types)):
        ax.bar(ctc_xranges[:-1]+x_offset, height=ctc_hists[ctc_types[i]], width=0.08,
               bottom=bottoms[i], label=ctc_types[i], color=colors[i])
    ax.set_xlim(ctc_xranges[0]-0.01, ctc_xranges[-1]+0.01)
    ax.set_ylim(0, 1)
    ax.set_xticks(ctc_xranges)
    ax.set_xticklabels(labels=[str(np.round(x,1)) for x in ctc_xranges],
                       fontdict={"fontsize":6})
    ax.legend(loc="upper right", fontsize=6)
    plt.ylabel("Fraction of contacts")
    plt.xlabel(r"Radii [$\AA$]")
    plt.title(r"Radii of closest sphere by contact residue type (3$\AA$) in 5-MDs or more")
    plt.savefig(os.path.join(outfolder, tag+"radii_by_restype_normalized_md5.png"), dpi=300, format="png")

def normalized_plot_HB_by_restype(data, hbond_outfolder,outfolder, tag="", frac_thr=0.7):
    bbone_atoms = {"N", "CA", "C", "O", "NH", "HA", "H"}
    aa_nonp = {"ALA", "GLY", "ILE", "LEU", "MET", "PHE", "PRO", "TRP", "VAL"}
    aa_polar = {"ASN", "CYS", "CYX", "GLN", "SER", "THR", "TYR"}
    aa_pos = {"ARG", "HID", "HIE", "HIP", "LYS"}
    aa_neg = {"GLU", "GLH", "ASP", "ASH"}
    
    def _get_hb_type(hbondfile):
        # Order of types:["bbone", "aa_nonp", "aa_polar", "aa_pos", "aa_neg"]
        restypes = [0, 0, 0, 0, 0]
        with open(hbondfile, "r") as hfile:
            [hfile.readline() for r in range(2)]
            for line in hfile:
                if line.startswith("#"):
                    break
                chunks = line.strip().split()
                if "Solvent" in chunks[0]:
                    res, atom = chunks[1].split("@")
                else:
                    res, atom = chunks[0].split("@")
                if atom in bbone_atoms:
                    restypes[0] += 1
                elif res[:3] in aa_nonp:
                    restypes[1] += 1
                elif res[:3] in aa_polar:
                    restypes[2] += 1
                elif res[:3] in aa_pos:
                    restypes[3] += 1
                elif res[:3] in aa_neg:
                    restypes[4] += 1
        return restypes
    
    mds = list(data.keys())
    mds.sort()
    for md in mds:
        for event in data[md]:
            hbondfile = "{}/{}/f{:0>5}_e{:0>4}_solv.txt".format(hbond_outfolder, md, event.frame+1, event.tid)
            event.restypes = _get_hb_type(hbondfile)
    
    sc_mds = defaultdict(int)
    for md in mds:
        _scs = set()
        for event in data[md]:
            if event.fraction >= frac_thr:
                _scs.add(event.supercluster)
        for s in _scs:
            sc_mds[s] += 1
    
    radii = defaultdict(list)
    restypes = ["bbone", "aa_nonp", "aa_polar", "aa_pos", "aa_neg"]
    colors = ["black", "tab:gray", "tab:olive", "tab:blue", "tab:red"]
    for md in mds:
        for event in data[md]:
            if (event.fraction >= frac_thr) and (sc_mds[event.supercluster] >=1):
                for i, rt in enumerate(restypes):
                    for j in range(event.restypes[i]):
                        radii[rt].append(event.radius)
    
    rt_hists = {}
    for rt in restypes:
        rt_hists[rt], hb_xranges = np.histogram(radii[rt], range=(0.9, 2.2), bins=10)
    x_offset = (hb_xranges[1] - hb_xranges[0]) / 2.0
    hist_sums = np.zeros(10, dtype=int)
    for rt in restypes:
        hist_sums += rt_hists[rt]
    for rt in restypes:
        rt_hists[rt] = np.nan_to_num(rt_hists[rt]/hist_sums, nan=0.0)
    
    bottoms = np.zeros((len(restypes), rt_hists["bbone"].shape[0]))
    for i in range(1, len(restypes)):
        for j in range(i, len(restypes)):
            bottoms[j,:] += rt_hists[restypes[i-1]]
    
    fig, ax = plt.subplots()
    for i in range(len(restypes)):
        ax.bar(hb_xranges[:-1]+x_offset, height=rt_hists[restypes[i]], width=0.09,
               bottom=bottoms[i], label=restypes[i], color=colors[i])
    ax.set_xlim(hb_xranges[0]-0.01, hb_xranges[-1]+0.01)
    ax.set_ylim(0, 1)
    ax.set_xticks(hb_xranges)
    ax.set_xticklabels(labels=[str(np.round(x,1)) for x in hb_xranges],
                       fontdict={"fontsize":6})
    ax.legend(loc="upper right", fontsize=6)
    plt.ylabel("Fraction of events")
    plt.xlabel(r"Radii [$\AA$]")
    plt.title("Radii of closest sphere to water transport by H-bond restype at least in 1-MDs")
    plt.savefig(os.path.join(outfolder, tag+"HB_by_restype_normalized_md1.png"), dpi=300, format="png")

def plot_radii(data, outfolder, tag="", frac_thr=0.7):
    mds = list(data.keys())
    mds.sort()
    radii_total = []
    
    sc_mds = defaultdict(int)
    for md in mds:
        _scs = set()
        for event in data[md]:
            if event.fraction >= frac_thr:
                _scs.add(event.supercluster)
        for s in _scs:
            sc_mds[s] += 1
    
    for md in mds:
        for event in data[md]:
            if (event.fraction >= frac_thr) and (sc_mds[event.supercluster] >= 1):
                radii_total.append(event.radius)
    
    fig, ax = plt.subplots()
    radii_hist, radii_xranges = np.histogram(radii_total, range=(0.7, 3), bins=30)
    x_offset = (radii_xranges[1] - radii_xranges[0]) / 2.0
    ax.bar(radii_xranges[:-1]+x_offset, height=radii_hist, width=0.07, color="tab:blue")
    ax.set_ylim(0, 80) # OPC : 25
    ax.set_xlim(radii_xranges[0]-0.01, radii_xranges[-1]+0.01)
    ax.set_xticks(radii_xranges)
    ax.set_xticklabels(labels=[str(np.round(x,1)) for x in radii_xranges],
                       fontdict={"fontsize":6})
    plt.ylabel("Number of events")
    plt.xlabel(r"Radii [$\AA$]")
    plt.title("Radii of closest sphere to water transport in 1-MDs or more")
    plt.savefig(os.path.join(outfolder, tag+"radii_md1.png"), dpi=300, format="png")

def plot_radii_by_SC(data, outfolder, tag="", frac_thr=0.7):
    mds = list(data.keys())
    mds.sort()
    radii = defaultdict(list)
    radii_hists = {}
    
    sc_mds = defaultdict(int)
    for md in mds:
        _scs = set()
        for event in data[md]:
            if event.fraction >= frac_thr:
                _scs.add(event.supercluster)
        for s in _scs:
            sc_mds[s] += 1
    
    for md in mds:
        for event in data[md]:
            if (event.fraction >= frac_thr) and (sc_mds[event.supercluster] >= 5):
                radii[event.supercluster].append(event.radius)
    
    scs = list(radii.keys())
    scs.sort()
    for sc in scs:
        radii_hists[sc], radii_xranges = np.histogram(radii[sc], range=(0.7, 3.7), bins=30)
    # radii_xranges = np.linspace(0.7, 3.7, 31)
    x_offset = (radii_xranges[1] - radii_xranges[0]) / 2.0
    bottoms = np.zeros((len(scs), radii_hists[scs[0]].shape[0]))
    for i in range(1, len(scs)):
        for j in range(i, len(scs)):
            bottoms[j,:] += radii_hists[scs[i-1]]
    
    fig, ax = plt.subplots()
    for i in range(len(scs)):
        ax.bar(radii_xranges[:-1]+x_offset, height=radii_hists[scs[i]], width=0.08,
               bottom=bottoms[i], label="SC_{}".format(scs[i]))
    ax.set_xlim(radii_xranges[0]-0.01, radii_xranges[-1]+0.01)
    ax.set_ylim(0, 1400) #dhaa:300, epx:2250, lipase:1400
    ax.set_xticks(radii_xranges)
    ax.set_xticklabels(labels=[str(np.round(x,1)) for x in radii_xranges],
                       fontdict={"fontsize":6})
    ax.legend(loc="upper right", fontsize=8)
    plt.ylabel("Number of events")
    plt.xlabel(r"Radii [$\AA$]")
    plt.title("Radii of closest sphere to water transport by SC in 5-MDs or more")
    plt.savefig(os.path.join(outfolder, tag+"radii_by_SC_md5.png"), dpi=300, format="png")

def plot_radii_by_HB(data,hbond_outfolder, outfolder, tag="", frac_thr=0.7):
    def _get_hbonds(hbondfile):
        hbonds = 0
        with open(hbondfile, "r") as hfile:
            for line in hfile:
                if line.startswith("#"):
                    continue
                hbonds = int(line.strip().split()[2])
        return hbonds
    
    mds = list(data.keys())
    mds.sort()
    for md in mds:
        for event in data[md]:
            # hbondfile = "{}/hbonds/{}/f{:0>5}_e{:0>4}.txt".format(outfolder, md, event.frame+1, event.tid)
            hbondfile = "{}/{}/f{:0>5}_e{:0>4}.txt".format(hbond_outfolder, md, event.frame + 1, event.tid)
            event.hbonds = _get_hbonds(hbondfile)
    
    sc_mds = defaultdict(int)
    for md in mds:
        _scs = set()
        for event in data[md]:
            if event.fraction >= frac_thr:
                _scs.add(event.supercluster)
        for s in _scs:
            sc_mds[s] += 1
    
    e_hb = defaultdict(list)
    for md in mds:
        for event in data[md]:
            if (event.fraction >= frac_thr) and (sc_mds[event.supercluster] >= 1):
                e_hb[event.hbonds].append(event.radius)
    
    ehb_hists = {}
    hbonds = list(e_hb.keys())
    hbonds.sort()
    for hb in hbonds:
        ehb_hists[hb], hb_xranges = np.histogram(e_hb[hb], range=(0.9, 3), bins=20)

    x_offset = (hb_xranges[1] - hb_xranges[0]) / 2.0
    
    bottoms = np.zeros((len(hbonds), ehb_hists[hbonds[0]].shape[0]))
    for i in range(1, len(hbonds)):
        for j in range(i, len(hbonds)):
            bottoms[j,:] += ehb_hists[hbonds[i-1]]
    
    fig, ax = plt.subplots()
    for i in range(len(hbonds)):
        ax.bar(hb_xranges[:-1]+x_offset, height=ehb_hists[hbonds[i]], width=0.09,
               bottom=bottoms[i], label="H-bonds_{}".format(hbonds[i]))
    ax.set_xlim(hb_xranges[0]-0.01, hb_xranges[-1]+0.01)
    ax.set_ylim(0,100) #dhaa:300, epx:2250, lipase:1400
    ax.set_xticks(hb_xranges)
    ax.set_xticklabels(labels=[str(np.round(x,1)) for x in hb_xranges],
                       fontdict={"fontsize":10})
    ax.legend(loc="upper right", fontsize=8)
    plt.ylabel("Number of events")
    plt.xlabel(r"Radii [$\AA$]")
    plt.title("Radii of closest sphere to water transport by H-bonds in 1-MDs or more")
    plt.savefig(os.path.join(outfolder, tag+"radii_by_HB_md1.png"), dpi=300, format="png")

def plot_radii_by_restype(data, outfolder, threshold=3.0, tag="", frac_thr=0.7):
    dhaa_res = {1:"ILE", 2:"GLY", 3:"THR", 4:"GLY", 5:"PHE", 6:"PRO", 7:"PHE", 8:"ASP", 9:"PRO", 10:"HID", 11:"TYR", 12:"VAL", 13:"GLU", 14:"VAL", 15:"LEU", 16:"GLY", 17:"GLU", 18:"ARG", 19:"MET", 20:"HID", 21:"TYR", 22:"VAL", 23:"ASP", 24:"VAL", 25:"GLY", 26:"PRO", 27:"ARG", 28:"ASP", 29:"GLY", 30:"THR", 31:"PRO", 32:"VAL", 33:"LEU", 34:"PHE", 35:"LEU", 36:"HID", 37:"GLY", 38:"ASN", 39:"PRO", 40:"THR", 41:"SER", 42:"SER", 43:"TYR", 44:"LEU", 45:"TRP", 46:"ARG", 47:"ASN", 48:"ILE", 49:"ILE", 50:"PRO", 51:"HIE", 52:"VAL", 53:"ALA", 54:"PRO", 55:"SER", 56:"HIE", 57:"ARG", 58:"CYS", 59:"ILE", 60:"ALA", 61:"PRO", 62:"ASP", 63:"LEU", 64:"ILE", 65:"GLY", 66:"MET", 67:"GLY", 68:"LYS", 69:"SER", 70:"ASP", 71:"LYS", 72:"PRO", 73:"ASP", 74:"LEU", 75:"ASP", 76:"TYR", 77:"PHE", 78:"PHE", 79:"ASP", 80:"ASP", 81:"HIE", 82:"VAL", 83:"ARG", 84:"TYR", 85:"LEU", 86:"ASP", 87:"ALA", 88:"PHE", 89:"ILE", 90:"GLU", 91:"ALA", 92:"LEU", 93:"GLY", 94:"LEU", 95:"GLU", 96:"GLU", 97:"VAL", 98:"VAL", 99:"LEU", 100:"VAL", 101:"ILE", 102:"HIE", 103:"ASP", 104:"TRP", 105:"GLY", 106:"SER", 107:"ALA", 108:"LEU", 109:"GLY", 110:"PHE", 111:"HIE", 112:"TRP", 113:"ALA", 114:"LYS", 115:"ARG", 116:"ASN", 117:"PRO", 118:"GLU", 119:"ARG", 120:"VAL", 121:"LYS", 122:"GLY", 123:"ILE", 124:"ALA", 125:"CYS", 126:"MET", 127:"GLU", 128:"PHE", 129:"ILE", 130:"ARG", 131:"PRO", 132:"ILE", 133:"PRO", 134:"THR", 135:"TRP", 136:"ASP", 137:"GLU", 138:"TRP", 139:"PRO", 140:"GLU", 141:"PHE", 142:"ALA", 143:"ARG", 144:"GLU", 145:"THR", 146:"PHE", 147:"GLN", 148:"ALA", 149:"PHE", 150:"ARG", 151:"THR", 152:"ALA", 153:"ASP", 154:"VAL", 155:"GLY", 156:"ARG", 157:"GLU", 158:"LEU", 159:"ILE", 160:"ILE", 161:"ASP", 162:"GLN", 163:"ASN", 164:"ALA", 165:"PHE", 166:"ILE", 167:"GLU", 168:"GLY", 169:"ALA", 170:"LEU", 171:"PRO", 172:"LYS", 173:"CYS", 174:"VAL", 175:"VAL", 176:"ARG", 177:"PRO", 178:"LEU", 179:"THR", 180:"GLU", 181:"VAL", 182:"GLU", 183:"MET", 184:"ASP", 185:"HIE", 186:"TYR", 187:"ARG", 188:"GLU", 189:"PRO", 190:"PHE", 191:"LEU", 192:"LYS", 193:"PRO", 194:"VAL", 195:"ASP", 196:"ARG", 197:"GLU", 198:"PRO", 199:"LEU", 200:"TRP", 201:"ARG", 202:"PHE", 203:"PRO", 204:"ASN", 205:"GLU", 206:"LEU", 207:"PRO", 208:"ILE", 209:"ALA", 210:"GLY", 211:"GLU", 212:"PRO", 213:"ALA", 214:"ASN", 215:"ILE", 216:"VAL", 217:"ALA", 218:"LEU", 219:"VAL", 220:"GLU", 221:"ALA", 222:"TYR", 223:"MET", 224:"ASN", 225:"TRP", 226:"LEU", 227:"HID", 228:"GLN", 229:"SER", 230:"PRO", 231:"VAL", 232:"PRO", 233:"LYS", 234:"LEU", 235:"LEU", 236:"PHE", 237:"TRP", 238:"GLY", 239:"THR", 240:"PRO", 241:"GLY", 242:"VAL", 243:"LEU", 244:"ILE", 245:"PRO", 246:"PRO", 247:"ALA", 248:"GLU", 249:"ALA", 250:"ALA", 251:"ARG", 252:"LEU", 253:"ALA", 254:"GLU", 255:"SER", 256:"LEU", 257:"PRO", 258:"ASN", 259:"CYS", 260:"LYS", 261:"THR", 262:"VAL", 263:"ASP", 264:"ILE", 265:"GLY", 266:"PRO", 267:"GLY", 268:"LEU", 269:"HIP", 270:"TYR", 271:"LEU", 272:"GLN", 273:"GLU", 274:"ASP", 275:"ASN", 276:"PRO", 277:"ASP", 278:"LEU", 279:"ILE", 280:"GLY", 281:"SER", 282:"GLU", 283:"ILE", 284:"ALA", 285:"ARG", 286:"TRP", 287:"LEU", 288:"PRO", 289:"ALA", 290:"LEU", 291:"HIE", 292:"HIE", 293:"HIP"}
    epx_res = {1:"LYS", 2:"ILE", 3:"GLU", 4:"HIE", 5:"LYS", 6:"MET", 7:"VAL", 8:"ALA", 9:"VAL", 10:"ASN", 11:"GLY", 12:"LEU", 13:"ASN", 14:"MET", 15:"HIP", 16:"LEU", 17:"ALA", 18:"GLU", 19:"LEU", 20:"GLY", 21:"GLU", 22:"GLY", 23:"PRO", 24:"THR", 25:"ILE", 26:"LEU", 27:"PHE", 28:"ILE", 29:"HID", 30:"GLY", 31:"PHE", 32:"PRO", 33:"GLU", 34:"LEU", 35:"TRP", 36:"TYR", 37:"SER", 38:"TRP", 39:"ARG", 40:"HIE", 41:"GLN", 42:"MET", 43:"VAL", 44:"TYR", 45:"LEU", 46:"ALA", 47:"GLU", 48:"ARG", 49:"GLY", 50:"TYR", 51:"ARG", 52:"ALA", 53:"VAL", 54:"ALA", 55:"PRO", 56:"ASP", 57:"LEU", 58:"ARG", 59:"GLY", 60:"TYR", 61:"GLY", 62:"ASP", 63:"THR", 64:"THR", 65:"GLY", 66:"ALA", 67:"PRO", 68:"LEU", 69:"ASN", 70:"ASP", 71:"PRO", 72:"SER", 73:"LYS", 74:"PHE", 75:"SER", 76:"ILE", 77:"LEU", 78:"HIE", 79:"LEU", 80:"VAL", 81:"GLY", 82:"ASP", 83:"VAL", 84:"VAL", 85:"ALA", 86:"LEU", 87:"LEU", 88:"GLU", 89:"ALA", 90:"ILE", 91:"ALA", 92:"PRO", 93:"ASN", 94:"GLU", 95:"GLU", 96:"LYS", 97:"VAL", 98:"PHE", 99:"VAL", 100:"VAL", 101:"ALA", 102:"HIE", 103:"ASP", 104:"TRP", 105:"GLY", 106:"ALA", 107:"LEU", 108:"ILE", 109:"ALA", 110:"TRP", 111:"HID", 112:"LEU", 113:"CYS", 114:"LEU", 115:"PHE", 116:"ARG", 117:"PRO", 118:"ASP", 119:"LYS", 120:"VAL", 121:"LYS", 122:"ALA", 123:"LEU", 124:"VAL", 125:"ASN", 126:"LEU", 127:"SER", 128:"VAL", 129:"HIE", 130:"PHE", 131:"SER", 132:"LYS", 133:"ARG", 134:"ASN", 135:"PRO", 136:"LYS", 137:"MET", 138:"ASN", 139:"VAL", 140:"VAL", 141:"GLU", 142:"GLY", 143:"LEU", 144:"LYS", 145:"ALA", 146:"ILE", 147:"TYR", 148:"GLY", 149:"GLU", 150:"ASP", 151:"HIE", 152:"TYR", 153:"ILE", 154:"SER", 155:"ARG", 156:"PHE", 157:"GLN", 158:"VAL", 159:"PRO", 160:"GLY", 161:"GLU", 162:"ILE", 163:"GLU", 164:"ALA", 165:"GLU", 166:"PHE", 167:"ALA", 168:"PRO", 169:"ILE", 170:"GLY", 171:"ALA", 172:"LYS", 173:"SER", 174:"VAL", 175:"LEU", 176:"LYS", 177:"LYS", 178:"ILE", 179:"LEU", 180:"THR", 181:"TYR", 182:"ARG", 183:"ASP", 184:"PRO", 185:"ALA", 186:"PRO", 187:"PHE", 188:"TYR", 189:"PHE", 190:"PRO", 191:"LYS", 192:"GLY", 193:"LYS", 194:"GLY", 195:"LEU", 196:"GLU", 197:"ALA", 198:"ILE", 199:"PRO", 200:"ASP", 201:"ALA", 202:"PRO", 203:"VAL", 204:"ALA", 205:"LEU", 206:"SER", 207:"SER", 208:"TRP", 209:"LEU", 210:"SER", 211:"GLU", 212:"GLU", 213:"GLU", 214:"LEU", 215:"ASP", 216:"TYR", 217:"TYR", 218:"ALA", 219:"ASN", 220:"LYS", 221:"PHE", 222:"GLU", 223:"GLN", 224:"THR", 225:"GLY", 226:"PHE", 227:"THR", 228:"GLY", 229:"ALA", 230:"VAL", 231:"ASN", 232:"TYR", 233:"TYR", 234:"ARG", 235:"ALA", 236:"LEU", 237:"PRO", 238:"ILE", 239:"ASN", 240:"TRP", 241:"GLU", 242:"LEU", 243:"THR", 244:"ALA", 245:"PRO", 246:"TRP", 247:"THR", 248:"GLY", 249:"ALA", 250:"GLN", 251:"VAL", 252:"LYS", 253:"VAL", 254:"PRO", 255:"THR", 256:"LYS", 257:"PHE", 258:"ILE", 259:"VAL", 260:"GLY", 261:"GLU", 262:"PHE", 263:"ASP", 264:"LEU", 265:"VAL", 266:"TYR", 267:"HIP", 268:"ILE", 269:"PRO", 270:"GLY", 271:"ALA", 272:"LYS", 273:"GLU", 274:"TYR", 275:"ILE", 276:"HIE", 277:"ASN", 278:"GLY", 279:"GLY", 280:"PHE", 281:"LYS", 282:"LYS", 283:"ASP", 284:"VAL", 285:"PRO", 286:"LEU", 287:"LEU", 288:"GLU", 289:"GLU", 290:"VAL", 291:"VAL", 292:"VAL", 293:"LEU", 294:"GLU", 295:"GLY", 296:"ALA", 297:"ALA", 298:"HIP", 299:"PHE", 300:"VAL", 301:"SER", 302:"GLN", 303:"GLU", 304:"ARG", 305:"PRO", 306:"HIE", 307:"GLU", 308:"ILE", 309:"SER", 310:"LYS", 311:"HIE", 312:"ILE", 313:"TYR", 314:"ASP", 315:"PHE", 316:"ILE", 317:"GLN", 318:"LYS", 319:"PHE"}
    lipase_res = {1:"ALA", 2:"PRO", 3:"THR", 4:"ALA", 5:"THR", 6:"LEU", 7:"ALA", 8:"ASN", 9:"GLY", 10:"ASP", 11:"THR", 12:"ILE", 13:"THR", 14:"GLY", 15:"LEU", 16:"ASN", 17:"ALA", 18:"ILE", 19:"ILE", 20:"ASN", 21:"GLU", 22:"ALA", 23:"PHE", 24:"LEU", 25:"GLY", 26:"ILE", 27:"PRO", 28:"PHE", 29:"ALA", 30:"GLU", 31:"PRO", 32:"PRO", 33:"VAL", 34:"GLY", 35:"ASN", 36:"LEU", 37:"ARG", 38:"PHE", 39:"LYS", 40:"ASP", 41:"PRO", 42:"VAL", 43:"PRO", 44:"TYR", 45:"SER", 46:"GLY", 47:"SER", 48:"LEU", 49:"ASP", 50:"GLY", 51:"GLN", 52:"LYS", 53:"PHE", 54:"THR", 55:"SER", 56:"TYR", 57:"GLY", 58:"PRO", 59:"SER", 60:"CYX", 61:"MET", 62:"GLN", 63:"GLN", 64:"ASN", 65:"PRO", 66:"GLU", 67:"GLY", 68:"THR", 69:"TYR", 70:"GLU", 71:"GLU", 72:"ASN", 73:"LEU", 74:"PRO", 75:"LYS", 76:"ALA", 77:"ALA", 78:"LEU", 79:"ASP", 80:"LEU", 81:"VAL", 82:"MET", 83:"GLN", 84:"SER", 85:"LYS", 86:"VAL", 87:"PHE", 88:"GLU", 89:"ALA", 90:"VAL", 91:"SER", 92:"PRO", 93:"SER", 94:"SER", 95:"GLU", 96:"ASP", 97:"CYX", 98:"LEU", 99:"THR", 100:"ILE", 101:"ASN", 102:"VAL", 103:"VAL", 104:"ARG", 105:"PRO", 106:"PRO", 107:"GLY", 108:"THR", 109:"LYS", 110:"ALA", 111:"GLY", 112:"ALA", 113:"ASN", 114:"LEU", 115:"PRO", 116:"VAL", 117:"MET", 118:"LEU", 119:"TRP", 120:"ILE", 121:"PHE", 122:"GLY", 123:"GLY", 124:"GLY", 125:"PHE", 126:"GLH", 127:"VAL", 128:"GLY", 129:"GLY", 130:"THR", 131:"SER", 132:"THR", 133:"PHE", 134:"PRO", 135:"PRO", 136:"ALA", 137:"GLN", 138:"MET", 139:"ILE", 140:"THR", 141:"LYS", 142:"SER", 143:"ILE", 144:"ALA", 145:"MET", 146:"GLY", 147:"LYS", 148:"PRO", 149:"ILE", 150:"ILE", 151:"HID", 152:"VAL", 153:"SER", 154:"VAL", 155:"ASN", 156:"TYR", 157:"ARG", 158:"VAL", 159:"SER", 160:"SER", 161:"TRP", 162:"GLY", 163:"PHE", 164:"LEU", 165:"ALA", 166:"GLY", 167:"ASP", 168:"GLU", 169:"ILE", 170:"LYS", 171:"ALA", 172:"GLU", 173:"GLY", 174:"SER", 175:"ALA", 176:"ASN", 177:"ALA", 178:"GLY", 179:"LEU", 180:"LYS", 181:"ASP", 182:"GLN", 183:"ARG", 184:"LEU", 185:"GLY", 186:"MET", 187:"GLN", 188:"TRP", 189:"VAL", 190:"ALA", 191:"ASP", 192:"ASN", 193:"ILE", 194:"ALA", 195:"ALA", 196:"PHE", 197:"GLY", 198:"GLY", 199:"ASP", 200:"PRO", 201:"THR", 202:"LYS", 203:"VAL", 204:"THR", 205:"ILE", 206:"PHE", 207:"GLY", 208:"GLU", 209:"SER", 210:"ALA", 211:"GLY", 212:"SER", 213:"MET", 214:"SER", 215:"VAL", 216:"MET", 217:"CYS", 218:"HIE", 219:"ILE", 220:"LEU", 221:"TRP", 222:"ASN", 223:"ASP", 224:"GLY", 225:"ASP", 226:"ASN", 227:"THR", 228:"TYR", 229:"LYS", 230:"GLY", 231:"LYS", 232:"PRO", 233:"LEU", 234:"PHE", 235:"ARG", 236:"ALA", 237:"GLY", 238:"ILE", 239:"MET", 240:"GLN", 241:"SER", 242:"GLY", 243:"ALA", 244:"MET", 245:"VAL", 246:"PRO", 247:"SER", 248:"ASP", 249:"ALA", 250:"VAL", 251:"ASP", 252:"GLY", 253:"ILE", 254:"TYR", 255:"GLY", 256:"ASN", 257:"GLU", 258:"ILE", 259:"PHE", 260:"ASP", 261:"LEU", 262:"LEU", 263:"ALA", 264:"SER", 265:"ASN", 266:"ALA", 267:"GLY", 268:"CYX", 269:"GLY", 270:"SER", 271:"ALA", 272:"SER", 273:"ASP", 274:"LYS", 275:"LEU", 276:"ALA", 277:"CYX", 278:"LEU", 279:"ARG", 280:"GLY", 281:"VAL", 282:"SER", 283:"SER", 284:"ASP", 285:"THR", 286:"LEU", 287:"GLH", 288:"ASP", 289:"ALA", 290:"THR", 291:"ASN", 292:"ASN", 293:"THR", 294:"PRO", 295:"GLY", 296:"PHE", 297:"LEU", 298:"ALA", 299:"TYR", 300:"SER", 301:"SER", 302:"LEU", 303:"ARG", 304:"LEU", 305:"SER", 306:"TYR", 307:"LEU", 308:"PRO", 309:"ARG", 310:"PRO", 311:"ASP", 312:"GLY", 313:"VAL", 314:"ASN", 315:"ILE", 316:"THR", 317:"ASP", 318:"ASP", 319:"MET", 320:"TYR", 321:"ALA", 322:"LEU", 323:"VAL", 324:"ARG", 325:"GLU", 326:"GLY", 327:"LYS", 328:"TYR", 329:"ALA", 330:"ASN", 331:"ILE", 332:"PRO", 333:"VAL", 334:"ILE", 335:"ILE", 336:"GLY", 337:"ASP", 338:"GLN", 339:"ASN", 340:"ASP", 341:"GLU", 342:"GLY", 343:"THR", 344:"PHE", 345:"PHE", 346:"GLY", 347:"THR", 348:"SER", 349:"SER", 350:"LEU", 351:"ASN", 352:"VAL", 353:"THR", 354:"THR", 355:"ASP", 356:"ALA", 357:"GLN", 358:"ALA", 359:"ARG", 360:"GLU", 361:"TYR", 362:"PHE", 363:"LYS", 364:"GLN", 365:"SER", 366:"PHE", 367:"VAL", 368:"HIP", 369:"ALA", 370:"SER", 371:"ASP", 372:"ALA", 373:"GLU", 374:"ILE", 375:"ASP", 376:"THR", 377:"LEU", 378:"MET", 379:"THR", 380:"ALA", 381:"TYR", 382:"PRO", 383:"GLY", 384:"ASP", 385:"ILE", 386:"THR", 387:"GLN", 388:"GLY", 389:"SER", 390:"PRO", 391:"PHE", 392:"ASP", 393:"THR", 394:"GLY", 395:"ILE", 396:"LEU", 397:"ASN", 398:"ALA", 399:"LEU", 400:"THR", 401:"PRO", 402:"GLN", 403:"PHE", 404:"LYS", 405:"ARG", 406:"ILE", 407:"SER", 408:"ALA", 409:"VAL", 410:"LEU", 411:"GLY", 412:"ASP", 413:"LEU", 414:"GLY", 415:"PHE", 416:"THR", 417:"LEU", 418:"ALA", 419:"ARG", 420:"ARG", 421:"TYR", 422:"PHE", 423:"LEU", 424:"ASN", 425:"HID", 426:"TYR", 427:"THR", 428:"GLY", 429:"GLY", 430:"THR", 431:"LYS", 432:"TYR", 433:"SER", 434:"PHE", 435:"LEU", 436:"SER", 437:"LYS", 438:"GLN", 439:"LEU", 440:"SER", 441:"GLY", 442:"LEU", 443:"PRO", 444:"VAL", 445:"LEU", 446:"GLY", 447:"THR", 448:"PHE", 449:"HIP", 450:"SER", 451:"ASN", 452:"ASP", 453:"ILE", 454:"VAL", 455:"PHE", 456:"GLN", 457:"ASH", 458:"TYR", 459:"LEU", 460:"LEU", 461:"GLY", 462:"SER", 463:"GLY", 464:"SER", 465:"LEU", 466:"ILE", 467:"TYR", 468:"ASN", 469:"ASN", 470:"ALA", 471:"PHE", 472:"ILE", 473:"ALA", 474:"PHE", 475:"ALA", 476:"THR", 477:"ASP", 478:"LEU", 479:"ASP", 480:"PRO", 481:"ASN", 482:"THR", 483:"ALA", 484:"GLY", 485:"LEU", 486:"LEU", 487:"VAL", 488:"LYS", 489:"TRP", 490:"PRO", 491:"GLU", 492:"TYR", 493:"THR", 494:"SER", 495:"SER", 496:"SER", 497:"GLN", 498:"SER", 499:"GLY", 500:"ASN", 501:"ASN", 502:"LEU", 503:"MET", 504:"MET", 505:"ILE", 506:"ASN", 507:"ALA", 508:"LEU", 509:"GLY", 510:"LEU", 511:"TYR", 512:"THR", 513:"GLY", 514:"LYS", 515:"ASP", 516:"ASN", 517:"PHE", 518:"ARG", 519:"THR", 520:"ALA", 521:"GLY", 522:"TYR", 523:"ASP", 524:"ALA", 525:"LEU", 526:"PHE", 527:"SER", 528:"ASN", 529:"PRO", 530:"PRO", 531:"SER", 532:"PHE", 533:"PHE", 534:"VAL"}
    hepx_wt_res = {1:"THR", 2:"SER", 3:"CYS", 4:"ASN", 5:"PRO", 6:"SER", 7:"ASP", 8:"MET", 9:"SER", 10:"HIE", 11:"GLY", 12:"TYR", 13:"VAL", 14:"THR", 15:"VAL", 16:"LYS", 17:"PRO", 18:"ARG", 19:"VAL", 20:"ARG", 21:"LEU", 22:"HIP", 23:"PHE", 24:"VAL", 25:"GLU", 26:"LEU", 27:"GLY", 28:"SER", 29:"GLY", 30:"PRO", 31:"ALA", 32:"VAL", 33:"CYS", 34:"LEU", 35:"CYS", 36:"HID", 37:"GLY", 38:"PHE", 39:"PRO", 40:"GLU", 41:"SER", 42:"TRP", 43:"TYR", 44:"SER", 45:"TRP", 46:"ARG", 47:"TYR", 48:"GLN", 49:"ILE", 50:"PRO", 51:"ALA", 52:"LEU", 53:"ALA", 54:"GLN", 55:"ALA", 56:"GLY", 57:"TYR", 58:"ARG", 59:"VAL", 60:"LEU", 61:"ALA", 62:"MET", 63:"ASP", 64:"MET", 65:"LYS", 66:"GLY", 67:"TYR", 68:"GLY", 69:"GLU", 70:"SER", 71:"SER", 72:"ALA", 73:"PRO", 74:"PRO", 75:"GLU", 76:"ILE", 77:"GLU", 78:"GLU", 79:"TYR", 80:"CYS", 81:"MET", 82:"GLU", 83:"VAL", 84:"LEU", 85:"CYS", 86:"LYS", 87:"GLU", 88:"MET", 89:"VAL", 90:"THR", 91:"PHE", 92:"LEU", 93:"ASP", 94:"LYS", 95:"LEU", 96:"GLY", 97:"LEU", 98:"SER", 99:"GLN", 100:"ALA", 101:"VAL", 102:"PHE", 103:"ILE", 104:"GLY", 105:"HIE", 106:"ASP", 107:"TRP", 108:"GLY", 109:"GLY", 110:"MET", 111:"LEU", 112:"VAL", 113:"TRP", 114:"TYR", 115:"MET", 116:"ALA", 117:"LEU", 118:"PHE", 119:"TYR", 120:"PRO", 121:"GLU", 122:"ARG", 123:"VAL", 124:"ARG", 125:"ALA", 126:"VAL", 127:"ALA", 128:"SER", 129:"LEU", 130:"ASN", 131:"THR", 132:"PRO", 133:"PHE", 134:"ILE", 135:"PRO", 136:"ALA", 137:"ASN", 138:"PRO", 139:"ASN", 140:"MET", 141:"SER", 142:"PRO", 143:"LEU", 144:"GLU", 145:"SER", 146:"ILE", 147:"LYS", 148:"ALA", 149:"ASN", 150:"PRO", 151:"VAL", 152:"PHE", 153:"ASP", 154:"TYR", 155:"GLN", 156:"LEU", 157:"TYR", 158:"PHE", 159:"GLN", 160:"GLU", 161:"PRO", 162:"GLY", 163:"VAL", 164:"ALA", 165:"GLU", 166:"ALA", 167:"GLU", 168:"LEU", 169:"GLU", 170:"GLN", 171:"ASN", 172:"LEU", 173:"SER", 174:"ARG", 175:"THR", 176:"PHE", 177:"LYS", 178:"SER", 179:"LEU", 180:"PHE", 181:"ARG", 182:"ALA", 183:"SER", 184:"ASP", 185:"GLU", 186:"SER", 187:"VAL", 188:"LEU", 189:"SER", 190:"MET", 191:"HIE", 192:"LYS", 193:"VAL", 194:"CYS", 195:"GLU", 196:"ALA", 197:"GLY", 198:"GLY", 199:"LEU", 200:"PHE", 201:"VAL", 202:"ASN", 203:"SER", 204:"PRO", 205:"GLU", 206:"GLU", 207:"PRO", 208:"SER", 209:"LEU", 210:"SER", 211:"ARG", 212:"MET", 213:"VAL", 214:"THR", 215:"GLU", 216:"GLU", 217:"GLU", 218:"ILE", 219:"GLN", 220:"PHE", 221:"TYR", 222:"VAL", 223:"GLN", 224:"GLN", 225:"PHE", 226:"LYS", 227:"LYS", 228:"SER", 229:"GLY", 230:"PHE", 231:"ARG", 232:"GLY", 233:"PRO", 234:"LEU", 235:"ASN", 236:"TRP", 237:"TYR", 238:"ARG", 239:"ASN", 240:"MET", 241:"GLU", 242:"ARG", 243:"ASN", 244:"TRP", 245:"LYS", 246:"TRP", 247:"ALA", 248:"CYS", 249:"LYS", 250:"SER", 251:"LEU", 252:"GLY", 253:"ARG", 254:"LYS", 255:"ILE", 256:"LEU", 257:"ILE", 258:"PRO", 259:"ALA", 260:"LEU", 261:"MET", 262:"VAL", 263:"THR", 264:"ALA", 265:"GLU", 266:"LYS", 267:"ASP", 268:"PHE", 269:"VAL", 270:"LEU", 271:"VAL", 272:"PRO", 273:"GLN", 274:"MET", 275:"SER", 276:"GLN", 277:"HIE", 278:"MET", 279:"GLU", 280:"ASP", 281:"TRP", 282:"ILE", 283:"PRO", 284:"HIE", 285:"LEU", 286:"LYS", 287:"ARG", 288:"GLY", 289:"HIE", 290:"ILE", 291:"GLU", 292:"ASP", 293:"CYS", 294:"GLY", 295:"HIP", 296:"TRP", 297:"THR", 298:"GLN", 299:"MET", 300:"ASP", 301:"LYS", 302:"PRO", 303:"THR", 304:"GLU", 305:"VAL", 306:"ASN", 307:"GLN", 308:"ILE", 309:"LEU", 310:"ILE", 311:"LYS", 312:"TRP", 313:"LEU", 314:"ASP", 315:"SER", 316:"ASP"}
    hepx_mut_res = {1:"THR", 2:"SER", 3:"CYS", 4:"ASN", 5:"PRO", 6:"SER", 7:"ASP", 8:"MET", 9:"SER", 10:"HIE", 11:"GLY", 12:"TYR", 13:"VAL", 14:"THR", 15:"VAL", 16:"LYS", 17:"PRO", 18:"ARG", 19:"VAL", 20:"ARG", 21:"LEU", 22:"HIP", 23:"PHE", 24:"VAL", 25:"GLU", 26:"LEU", 27:"GLY", 28:"SER", 29:"GLY", 30:"PRO", 31:"ALA", 32:"VAL", 33:"CYS", 34:"LEU", 35:"CYS", 36:"HID", 37:"GLY", 38:"PHE", 39:"PRO", 40:"GLU", 41:"SER", 42:"TRP", 43:"TYR", 44:"SER", 45:"TRP", 46:"ARG", 47:"TYR", 48:"GLN", 49:"ILE", 50:"PRO", 51:"ALA", 52:"LEU", 53:"ALA", 54:"GLN", 55:"ALA", 56:"GLY", 57:"TYR", 58:"ARG", 59:"VAL", 60:"LEU", 61:"ALA", 62:"MET", 63:"ASP", 64:"MET", 65:"LYS", 66:"GLY", 67:"TYR", 68:"GLY", 69:"GLU", 70:"SER", 71:"SER", 72:"ALA", 73:"PRO", 74:"PRO", 75:"GLU", 76:"ILE", 77:"GLU", 78:"GLU", 79:"TYR", 80:"CYS", 81:"MET", 82:"GLU", 83:"VAL", 84:"LEU", 85:"CYS", 86:"LYS", 87:"GLU", 88:"MET", 89:"VAL", 90:"THR", 91:"PHE", 92:"LEU", 93:"ASP", 94:"LYS", 95:"LEU", 96:"GLY", 97:"LEU", 98:"SER", 99:"GLN", 100:"ALA", 101:"VAL", 102:"PHE", 103:"ILE", 104:"GLY", 105:"HIE", 106:"ASP", 107:"TRP", 108:"GLY", 109:"GLY", 110:"MET", 111:"LEU", 112:"VAL", 113:"TRP", 114:"TYR", 115:"MET", 116:"ALA", 117:"LEU", 118:"PHE", 119:"TYR", 120:"PRO", 121:"GLU", 122:"ARG", 123:"VAL", 124:"ARG", 125:"ALA", 126:"VAL", 127:"ALA", 128:"SER", 129:"LEU", 130:"ASN", 131:"THR", 132:"PRO", 133:"PHE", 134:"ILE", 135:"PRO", 136:"ALA", 137:"ASN", 138:"PRO", 139:"ASN", 140:"MET", 141:"SER", 142:"PRO", 143:"LEU", 144:"GLU", 145:"SER", 146:"ILE", 147:"LYS", 148:"ALA", 149:"ASN", 150:"PRO", 151:"VAL", 152:"PHE", 153:"ASP", 154:"TYR", 155:"GLN", 156:"LEU", 157:"TYR", 158:"PHE", 159:"GLN", 160:"GLU", 161:"PRO", 162:"GLY", 163:"VAL", 164:"ALA", 165:"GLU", 166:"ALA", 167:"GLU", 168:"LEU", 169:"GLU", 170:"GLN", 171:"ASN", 172:"LEU", 173:"SER", 174:"ARG", 175:"THR", 176:"PHE", 177:"LYS", 178:"SER", 179:"LEU", 180:"PHE", 181:"ARG", 182:"ALA", 183:"SER", 184:"ASP", 185:"GLU", 186:"SER", 187:"VAL", 188:"LEU", 189:"SER", 190:"MET", 191:"HIE", 192:"LYS", 193:"VAL", 194:"CYS", 195:"GLU", 196:"ALA", 197:"GLY", 198:"GLY", 199:"LEU", 200:"PHE", 201:"VAL", 202:"ASN", 203:"SER", 204:"PRO", 205:"GLU", 206:"GLU", 207:"PRO", 208:"SER", 209:"LEU", 210:"SER", 211:"ARG", 212:"MET", 213:"VAL", 214:"THR", 215:"GLU", 216:"GLU", 217:"GLU", 218:"ILE", 219:"GLN", 220:"PHE", 221:"TYR", 222:"VAL", 223:"GLN", 224:"GLN", 225:"PHE", 226:"LYS", 227:"LYS", 228:"SER", 229:"GLY", 230:"PHE", 231:"ARG", 232:"GLY", 233:"PRO", 234:"LEU", 235:"ASN", 236:"TRP", 237:"TYR", 238:"ARG", 239:"ASN", 240:"MET", 241:"GLY", 242:"ARG", 243:"ASN", 244:"TRP", 245:"LYS", 246:"TRP", 247:"ALA", 248:"CYS", 249:"LYS", 250:"SER", 251:"LEU", 252:"GLY", 253:"ARG", 254:"LYS", 255:"ILE", 256:"LEU", 257:"ILE", 258:"PRO", 259:"ALA", 260:"LEU", 261:"MET", 262:"VAL", 263:"THR", 264:"ALA", 265:"GLU", 266:"LYS", 267:"ASP", 268:"PHE", 269:"VAL", 270:"LEU", 271:"VAL", 272:"PRO", 273:"GLN", 274:"MET", 275:"SER", 276:"GLN", 277:"HIE", 278:"MET", 279:"GLU", 280:"ASP", 281:"TRP", 282:"ILE", 283:"PRO", 284:"HIE", 285:"LEU", 286:"LYS", 287:"ARG", 288:"GLY", 289:"HIE", 290:"ILE", 291:"GLU", 292:"ASP", 293:"CYS", 294:"GLY", 295:"HIP", 296:"TRP", 297:"THR", 298:"GLN", 299:"MET", 300:"ASP", 301:"LYS", 302:"PRO", 303:"THR", 304:"GLU", 305:"VAL", 306:"ASN", 307:"GLN", 308:"ILE", 309:"LEU", 310:"ILE", 311:"LYS", 312:"TRP", 313:"LEU", 314:"ASP", 315:"SER", 316:"ASP"}
    aa_neg = {"GLU", "GLH", "ASP", "ASH"}
    aa_pos = {"ARG", "HID", "HIE", "HIP", "LYS"}
    aa_polar = {"ASN", "CYS", "CYX", "GLN", "SER", "THR", "TYR"}
    aa_nonp = {"ALA", "GLY", "ILE", "LEU", "MET", "PHE", "PRO", "TRP", "VAL"}
    bbone_atoms = {"N", "CA", "C", "O", "NH", "HA", "H"}
    
    mds = list(data.keys())
    mds.sort()
    
    sc_mds = defaultdict(int)
    for md in mds:
        _scs = set()
        for event in data[md]:
            if event.fraction >= frac_thr:
                _scs.add(event.supercluster)
        for s in _scs:
            sc_mds[s] += 1
    
    residues = defaultdict(int)
    radii_by_type = defaultdict(list)
    prot = lipase_res
    for md in mds:
        for event in data[md]:
            _resids = set()
            contactfile = "{}/ncontacts/{}/f{:0>5}_e{:0>4}.txt".format(outfolder, md, event.frame+1, event.tid)
            if sc_mds[event.supercluster] >= 5:
                with open(contactfile, "r") as fin:
                    [fin.readline() for i in range(3)]
                    for line in fin:
                        chunks = line.strip().split()
                        dist = float(chunks[4])
                        _res = chunks[1][1:].split(":")[-1]
                        _res, _atom = _res.split("@")
                        if dist <= threshold:
                            r = int(_res)
                            if r not in _resids:
                                _resids.add(r)
                                if _atom in bbone_atoms:
                                    radii_by_type["bbone"].append(event.radius)
                                elif prot[r] in aa_neg:
                                    radii_by_type["aa_neg"].append(event.radius)
                                elif prot[r] in aa_pos:
                                    radii_by_type["aa_pos"].append(event.radius)
                                elif prot[r] in aa_polar:
                                    radii_by_type["aa_polar"].append(event.radius)
                                elif prot[r] in aa_nonp:
                                    radii_by_type["aa_nonp"].append(event.radius)
                                else:
                                    print("Residue not in groups:", prot[r])
    
    ctc_hists = {}
    ctc_types = ["bbone", "aa_nonp", "aa_polar", "aa_pos", "aa_neg"]
    colors = ["black", "tab:gray", "tab:olive", "tab:blue", "tab:red"]
    for ctc_type in ctc_types:
        ctc_hists[ctc_type], ctc_xranges = np.histogram(radii_by_type[ctc_type], range=(0.7, 3.7), bins=30)
    x_offset = (ctc_xranges[1] - ctc_xranges[0]) / 2.0
    
    bottoms = np.zeros((len(ctc_types), ctc_hists[ctc_types[0]].shape[0]))
    for i in range(1, len(ctc_types)):
        for j in range(i, len(ctc_types)):
            bottoms[j,:] += ctc_hists[ctc_types[i-1]]
    
    fg, ax = plt.subplots()
    for i in range(len(ctc_types)):
        ax.bar(ctc_xranges[:-1]+x_offset, height=ctc_hists[ctc_types[i]], width=0.08,
               bottom=bottoms[i], label=ctc_types[i], color=colors[i])
    ax.set_xlim(ctc_xranges[0]-0.01, ctc_xranges[-1]+0.01)
    # ax.set_ylim(0, 1400) #dhaa:300, epx:2250, lipase:1400
    ax.set_xticks(ctc_xranges)
    ax.set_xticklabels(labels=[str(np.round(x,1)) for x in ctc_xranges],
                       fontdict={"fontsize":6})
    ax.legend(loc="upper right", fontsize=8)
    plt.ylabel("Number of contacts")
    plt.xlabel(r"Radii [$\AA$]")
    plt.title("Radii of closest sphere to water transport by residue type in 5-MDs or more")
    plt.savefig(os.path.join(outfolder, tag+"radii_by_restype_md5.png"), dpi=300, format="png")

if __name__ == "__main__":
    # python3 10_new_histograms.py dhaa/events_data.dat dhaa dhaa
    # if (len(argv) != 4) and (len(argv) != 5):
    #     print("Usage:python3 05_analyze_ncontacts_hbonds.py event_datafile hbond_results_dir results_folder tag_for_figures(optional)")
    #     print("Usage:python3 05_analyze_ncontacts_hbonds.py dhaa/events_data.dat dhaa dhaa_figs")
    #     quit(1)
    # elif len(argv) == 5:
    #     tag = argv[4] + "_"
    # else:
    #     tag = ""
    # with open(argv[1], "rb") as fin:
    #     database = pickle.load(fin)
    # hbond_outfolder =argv[2]
    dat_file = "/data/aravindramt/dean/tt/hbond/databases/database_1_O.dat"
    with open(dat_file, "rb") as fin:
        database = pickle.load(fin)
    hbond_outfolder= "/data/aravindramt/dean/tt/hbond/hbonds/hb_1_O"
    plot_folder = "/data/aravindramt/dean/tt/hbond/plots/1_OPC"
    # os.mkdir(plot_folder)
    # tag = "1_4_OPC_"
    # normalized_plot_radii_by_SC(database, argv[3], tag=tag)

    # normalized_plot_radii_by_HB(database,hbond_outfolder, argv[3], tag=tag)


    # normalized_plot_radii_by_restype(database, argv[2], tag=tag)

    # normalized_plot_HB_by_restype(database,hbond_outfolder, argv[3], tag=tag)


    # plot_radii(database, argv[3], tag=tag)

    # plot_radii_by_SC(database, argv[3], tag=tag)
    # plot_radii_by_HB(database,hbond_outfolder, argv[3], tag=tag)

    # plot_radii_by_restype(database, argv[2], tag=tag)
    # plot_radii(database, plot_folder)
    plot_radii_by_HB(database, hbond_outfolder, plot_folder)
    # normalized_plot_HB_by_restype(database, hbond_outfolder, plot_folder)
    # normalized_plot_radii_by_HB(database, hbond_outfolder, plot_folder)
