import numpy as np
import os
import time
from vocustools.auxilary_routines import getpar, getstr, get_inst_id, get_tr, add_tr, formula, d, drucksieb, m2t, str2vec, t09, t09str
from vocustools.data_processing import *
from vocustools import *
from vocustools.auxilary_routines import *

def uni_mass_list(files, names, dest_folder):
    lib = masslib(extended=True)
    
    num_files = len(files)
    time_process = 0
    time_per_file = 15
    new_file = 0
    print(ReadFloat(dest_folder + 'Lastcheck.fdt'))
    last_check = max(ReadFloat(dest_folder + 'Lastcheck.fdt'))
    date_string = '01/01/1970, 00:00:00'
    
    now = max([t09(date_string) + time.time() / (24.0 * 3600.0)])
    
    exsts1 = os.path.exists(dest_folder + 'binScale.fdt')
    exsts2 = os.path.exists(dest_folder + 'resolution.fdt')
    exsts3 = os.path.exists(dest_folder + 'PeakPar.fdt')
    if isinstance(last_check, int):
        last_check = [last_check]
    if now - last_check[0] > 0.5 or exsts1 + exsts2 + exsts3 != 3:
        makefloat(dest_folder + 'Lastcheck.fdt', [now])
        if exsts1 + exsts2 + exsts3 != 3:
            new_file = 1
        
        for i in range(num_files):
            exists = os.path.exists(dest_folder + 'FileInfo/IonList/' + names[i] + 'FL2.fdt')
            if not exists:
                new_file = 1
                runtime_strt = time.time()
                eng_data = get_eng_data(files[i])
                
                if eng_data['duration'] > 0:
                    if 'TOF1000' in eng_data['instrument']:
                        min_res = 500
                    elif 'VOCUS' in eng_data['instrument']:
                        min_res = 5000
                    else:
                        min_res = 1800
                    
                    jj = 0
                    duration = eng_data['duration']
                    sumspectr = eng_data['sumspec']
                    next_start_time = eng_data['starttime'] + duration / (3600 * 24)
                    sig = max(np.convolve(sumspectr, np.ones(3)/3, mode='same'))
                    
                    while sig < getpar('Desired_Min_Signal'):
                        if i + jj + 1 >= num_files:
                            break
                        next_data = get_eng_data(files[i + jj + 1])
                        if next_data['duration'] == -9999:
                            break
                        if abs(next_data['starttime'] - next_start_time) * 24 * 60 > getpar('Max_Time_Gap'):
                            break
                        jj += 1
                        duration += next_data['duration']
                        sumspectr += next_data['sumspec']
                        sig += max(np.convolve(next_data['sumspec'], np.ones(3)/3, mode='same'))
                        next_start_time = next_data['starttime'] + next_data['duration'] / (3600 * 24)
                    
                    runtime_open_file = time.time() - runtime_strt
                    
                    if duration > 0:
                        good_f = 0
                        fls = [names[i]]
                        if jj > 0.5:
                            fls.extend(names[i+1:i+jj+1])
                        
                        print(files[i])
                        file_par = PS1(sumspectr, eng_data['SampInt'], eng_data['duration'], eng_data['cycles'], eng_data['extractions'], dest_folder, fls, eng_data['instrument'])
                        a_init = file_par['a']
                        t0_init = file_par['t0']
                        
                        if getpar('DriftCorrection') == 1:
                            sumspec = ReadFloat(dest_folder + 'FileInfo/Sumspec' + names[i] + '.fdt')
                            if max(sumspec) == -9999:
                                haha = getSplit(files[i])
                                split = haha['split']
                                misst = LoadCycles(files[i], haha['locatie'], haha['start2'][:, 0], haha['count'][:, 0])
                                de_dri = DeDrift(misst, a_init, t0_init, eng_data['SampInt'], -9999, -9999, -9999, -9999, -9999, -9999, split, 1)
                                sumspec = de_dri['sumspec']
                                if split > 1.5:
                                    for ppp in range(1, split):
                                        misst = LoadCycles(files[i], haha['locatie'], haha['start2'][:, ppp], haha['count'][:, ppp])
                                        de_dri = DeDrift(misst, a_init, t0_init, eng_data['SampInt'], de_dri['rest'], de_dri['maxis'], de_dri['t1start'], de_dri['t2start'], de_dri['m1'], de_dri['m2'], split, ppp + 1)
                                        sumspec += de_dri['sumspec']
                                makefloat(dest_folder + 'FileInfo/Sumspec' + names[i] + '.fdt', sumspec)
                            file_par = PS1(sumspec, eng_data['SampInt'], eng_data['duration'], eng_data['cycles'], eng_data['extractions'], dest_folder, fls, eng_data['instrument'])
                        
                        if max(file_par['counts']) != -9999:
                            file_par['counts'][file_par['counts'] > 0] /= duration
                        if max(file_par['res']) > min_res:
                            good_f = 1
                        
                        if max(file_par['a']) != -9999:
                            time_process += time.time() - runtime_strt
                            time_per_file.append(time.time() - runtime_strt)
                            detected_peaks = len(file_par['masslist'])
                            
                            # Log processing results
                            log = [
                                fls[0],
                                f"StrtTime [T09]: {eng_data['starttime']:.3f}",
                                f"Duration [min]: {duration/60:.2f}",
                                f"MaxSignal: {max(sumspectr):.2e}",
                                f"MaxSignal/s: {max(sumspectr)/duration:.0f}",
                                f"Maxmass [Da]: {file_par['maxmass']}",
                                f"a: {file_par['a']}",
                                f"t0: {file_par['t0']}",
                                f"ex: {file_par['ex']}",
                                f"Resolution (FWHM): {file_par['resolution']:5d}",
                                f"Detected peaks: {detected_peaks}",
                                "------------",
                                "",
                                "Test: mass, [cps], [FWHM]"
                            ]
                            
                            # Add test data
                            for ij in range(10):
                                log.append(f"{file_par['massnames'][ij]} {file_par['counts'][ij]:10.0f} {file_par['res'][ij]:10.0f}")
                            
                            log.extend([
                                f"Total processing time: {time.time() - runtime_strt:.1f}s",
                                f"Files remaining: {num_files - i - 1:4d}",
                                f"Time remaining: {np.quantile(time_per_file, 0.8) * (num_files - i - 1):6.0f}s",
                                "------------",
                                ""
                            ])
                            
                            # Truncate log if too long
                            if len(log) > 100:
                                log = log[:100]
                        
                        else:
                            log = [fls[0], "Crude calibration FAILED", "------------", ""]
                        
                        # Save processing results
                        for kk in range(jj + 1):
                            if good_f == 1:
                                makefloat(dest_folder + 'FileInfo/IonList/' + names[i+kk] + 'FL2.fdt', file_par['masslist'])
                                makefloat(dest_folder + 'FileInfo/Time2MassPAR/' + names[i+kk] + 'PAR2.fdt', [file_par['a'], file_par['t0'], file_par['ex'], file_par['maxmass'], eng_data['SampInt'], file_par['a3'], file_par['t03'], file_par['ex3'], file_par['resolution'], file_par['mode'].index('H3O+ NO+ O2+'), a_init, t0_init])
                                makefloat(dest_folder + 'FileInfo/PeakShape/' + names[i+kk] + 'PeakShape.fdt', file_par['PeakShape'])
                                makefloat(dest_folder + 'FileInfo/Baseline/' + names[i+kk] + 'Baseline.fdt', file_par['baseline'])
                            else:
                                makefloat(dest_folder + 'FileInfo/IonList/' + names[i+kk] + 'FL2.fdt', [-9999])
                                makefloat(dest_folder + 'FileInfo/Time2MassPAR/' + names[i+kk] + 'PAR2.fdt', [-9999])
                                makefloat(dest_folder + 'FileInfo/PeakShape/' + names[i+kk] + 'PeakShape.fdt', [-9999])
                                makefloat(dest_folder + 'FileInfo/Baseline/' + names[i+kk] + 'Baseline.fdt', [-9999])
                        i += jj
    
    runtime_strt = time.time()
    print("new_file", new_file, num_files )
    if new_file == 1:
        masslist = []
        resolution = []
        for i in range(num_files):
            if os.path.exists(dest_folder + 'FileInfo/IonList/' + names[i] + 'FL2.fdt'):
                FL2 = ReadFloat(dest_folder + 'FileInfo/IonList/' + names[i] + 'FL2.fdt')
                PAR2 = ReadFloat(dest_folder + 'FileInfo/Time2MassPAR/' + names[i] + 'PAR2.fdt')
                print("FL2", FL2,  PAR2)
                if PAR2[0] > 0:
                    resolution.append(PAR2[8])
                if FL2[0] > 0:
                    masslist.extend(FL2[:, 0])
            else: print("File not found: " + dest_folder + 'FileInfo/IonList/' + names[i] + 'FL2.fdt')
        
        binwidth = 32 if np.median(resolution) < 2500 else 8
        binscale = ppm_bin(masslist, binwidth)
        x = binscale[:, 0]
        y = binscale[:, 2]
        peakdata = PS2(x, y)
        masslist = peakdata['data'][:, 1]
        intlist = PeakTable(masslist, np.median(resolution))
        peakdat = peakdata['data']
        
        filter = intlist[:, 2] > 0.007 * 3000.0 / np.median(resolution) + intlist[:, 1]
        masslist = masslist[filter]
        intlist = intlist[filter]
        peakdat = peakdat[filter]
        
        makefloat(dest_folder + 'binScale.fdt', binscale)
        makefloat(dest_folder + 'resolution.fdt', resolution)
        MakeCsv(dest_folder + 'UnifiedMasslist.csv', intlist.T)
        MakeCsv(dest_folder + 'UnifiedMasslistNames.csv', ['mass[Da]', 'StartInt', 'EndInt', 'sigma', 'StartCut [fraction of sigma]', 'EndCut [fraction of sigma]'])
        makefloat(dest_folder + 'PeakPar.fdt', peakdat)
        MakeCsv(dest_folder + 'PeakPar.csv', peakdat.T)
        MakeCsv(dest_folder + 'PeakParNames.csv', peakdata['names'])
    else:
        binscale = ReadFloat(dest_folder + 'binScale.fdt')
        resolution = ReadFloat(dest_folder + 'resolution.fdt')
        x = binscale[:, 0]
        y = binscale[:, 2]
        peakdat = ReadFloat(dest_folder + 'PeakPar.fdt')
    
    xpeaks_der = peakdat[:, 0]
    xpeaks = peakdat[:, 1]
    sigma = peakdat[:, 2]
    ymax = peakdat[:, 3]
    
    return {
        'PeakCountsPERppmBINS': y,
        'PpmBINscale': x,
        'PeakLoc': xpeaks,
        'PeakLocDer': xpeaks_der,
        'PeakWidth': sigma,
        'PeakHeight': ymax
    }

def masslib(extended=False, sulfate=False, cal=False):
    # Define atomic masses
    C12, C13 = 12, 13.003355
    H1 = 1.007825
    O16, O17, O18 = 15.994915, 16.999131, 17.99916
    N14, N15 = 14.003074, 15.000108
    H_pos, e_neg = 1.007276467, -0.000548533
    Cl35, Cl37 = 34.968852, 36.965903
    F = 18.998403
    S32O4, S34O4 = 31.97207 + 4*O16, 33.967866 + 4*O16
    Si28, Si29, Si30 = 27.976927, 28.976495, 29.97377
    I = 126.904473

    # Set parameters based on keywords
    ooo = 16 if extended else 5
    nnn = 2 if extended else 0
    Clll = 0
    sss = 1 if sulfate else 0
    SRI = 1 if cal else 0

    if cal:
        ooo, nnn, Clll, sss = 5, 0, 0, 0

    # Initialize library
    lib = []

    # Add primary/inorganic ions
    primary_ions = [
        [0, 0, 2, 0, 0, 1, 0, 0, 1, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0],  # m21
        [0, 0, 4, 1, 0, 1, 0, 0, 1, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0],  # m39
        [0, 0, 0, 2, 0, 0, 0, 0, 0, 1, 0, 0, 0, 0, 0, 0, 0, 0, 0],  # O2+
        [0, 0, 0, 1, 0, 1, 0, 0, 0, 1, 0, 0, 0, 0, 0, 0, 0, 0, 0],  # 16O18O+
        [0, 0, 0, 1, 0, 0, 1, 0, 0, 1, 0, 0, 0, 0, 0, 0, 0, 0, 0],  # NO+
        [0, 0, 0, 1, 0, 0, 0, 1, 0, 1, 0, 0, 0, 0, 0, 0, 0, 0, 0],  # 15NO+
        [0, 0, 0, 0, 0, 0, 2, 0, 1, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0],  # N2H+
        [0, 0, 3, 0, 0, 0, 1, 0, 1, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0],  # NH4+
        [0, 0, 6, 0, 0, 0, 2, 0, 1, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0],  # NH3NH4+
        [0, 0, 1, 3, 0, 0, 1, 0, 1, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0],  # HNO3H+
        [0, 0, 0, 2, 0, 0, 1, 0, 0, 1, 0, 0, 0, 0, 0, 0, 0, 0, 0],  # NO2+
        [0, 0, 0, 2, 0, 0, 0, 1, 0, 1, 0, 0, 0, 0, 0, 0, 0, 0, 0],  # 15NO2+
        [0, 0, 0, 1, 0, 1, 1, 0, 0, 1, 0, 0, 0, 0, 0, 0, 0, 0, 0],  # NO18O+
    ]

    for n in primary_ions:
        mass = sum([n[i] * x for i, x in enumerate([C12, C13, H1, O16, O17, O18, N14, N15, H_pos, e_neg, Cl35, Cl37, F, S32O4, S34O4, Si28, Si29, Si30, I])])
        lib.append([mass] + n)

    # Add CCl3+ and its isotopes
    for c, cl35, cl37 in [(1, 3, 0), (1, 2, 1), (0, 3, 0), (0, 2, 1), (1, 1, 2)]:
        n = [c, 1-c, 0, 0, 0, 0, 0, 0, 0, 1, cl35, cl37, 0, 0, 0, 0, 0, 0, 0]
        mass = sum([n[i] * x for i, x in enumerate([C12, C13, H1, O16, O17, O18, N14, N15, H_pos, e_neg, Cl35, Cl37, F, S32O4, S34O4, Si28, Si29, Si30, I])])
        lib.append([mass] + n)

    # Add Iodobenzene and Di-Iodobenzene
    iodo_compounds = [
        [6, 0, 5, 0, 0, 0, 0, 0, 1, 0, 0, 0, 0, 0, 0, 0, 0, 0, 1],
        [5, 1, 5, 0, 0, 0, 0, 0, 1, 0, 0, 0, 0, 0, 0, 0, 0, 0, 1],
        [6, 0, 4, 0, 0, 0, 0, 0, 1, 0, 0, 0, 0, 0, 0, 0, 0, 0, 2],
        [5, 1, 4, 0, 0, 0, 0, 0, 1, 0, 0, 0, 0, 0, 0, 0, 0, 0, 2],
        [6, 0, 4, 0, 0, 0, 0, 0, 1, 0, 0, 0, 0, 0, 0, 0, 0, 0, 1],
        [5, 1, 4, 0, 0, 0, 0, 0, 1, 0, 0, 0, 0, 0, 0, 0, 0, 0, 1],
    ]

    for n in iodo_compounds:
        mass = sum([n[i] * x for i, x in enumerate([C12, C13, H1, O16, O17, O18, N14, N15, H_pos, e_neg, Cl35, Cl37, F, S32O4, S34O4, Si28, Si29, Si30, I])])
        lib.append([mass] + n)

    # Add 3F-benzene and 3Cl-benzene
    halogen_benzenes = [
        [6, 0, 3, 0, 0, 0, 0, 0, 1, 0, 0, 0, 3, 0, 0, 0, 0, 0, 0],
        [5, 1, 3, 0, 0, 0, 0, 0, 1, 0, 0, 0, 3, 0, 0, 0, 0, 0, 0],
        [6, 0, 3, 0, 0, 0, 0, 0, 1, 0, 3, 0, 0, 0, 0, 0, 0, 0, 0],
        [6, 0, 3, 0, 0, 0, 0, 0, 1, 0, 2, 1, 0, 0, 0, 0, 0, 0, 0],
        [5, 1, 3, 0, 0, 0, 0, 0, 1, 0, 3, 0, 0, 0, 0, 0, 0, 0, 0],
        [5, 1, 3, 0, 0, 0, 0, 0, 1, 0, 2, 1, 0, 0, 0, 0, 0, 0, 0],
        [6, 0, 3, 0, 0, 0, 0, 0, 1, 0, 1, 2, 0, 0, 0, 0, 0, 0, 0],
    ]

    for n in halogen_benzenes:
        mass = sum([n[i] * x for i, x in enumerate([C12, C13, H1, O16, O17, O18, N14, N15, H_pos, e_neg, Cl35, Cl37, F, S32O4, S34O4, Si28, Si29, Si30, I])])
        lib.append([mass] + n)

    # Generate D and L series compounds
    for ss in range(4):
        for series in ['D', 'L']:
            base_c = 6 if series == 'D' else 8
            base_h = 18 if series == 'D' else 24
            base_o = 3 if series == 'D' else 2

            for isotope in range(3):
                n = [base_c + ss*2, 0, base_h + ss*6, base_o + ss, 0, 0, 0, 0, 1, 0, 0, 0, 0, 0, 0, 3-isotope+ss, isotope, 0, 0]
                mass = sum([n[i] * x for i, x in enumerate([C12, C13, H1, O16, O17, O18, N14, N15, H_pos, e_neg, Cl35, Cl37, F, S32O4, S34O4, Si28, Si29, Si30, I])])
                lib.append([mass] + n)

            # Fragment -CH4
            for isotope in range(3):
                n = [base_c - 1 + ss*2, 0, base_h - 4 + ss*6, base_o + ss, 0, 0, 0, 0, 1, 0, 0, 0, 0, 0, 0, 3-isotope+ss, isotope, 0, 0]
                mass = sum([n[i] * x for i, x in enumerate([C12, C13, H1, O16, O17, O18, N14, N15, H_pos, e_neg, Cl35, Cl37, F, S32O4, S34O4, Si28, Si29, Si30, I])])
                lib.append([mass] + n)

    # Generate CHNOCl compounds
    for c in range(1, 41):
        for h in range(c-9, c+2):
            for o in range(ooo+1):
                for nn in range(nnn+1):
                    for Cll in range(Clll+1):
                        for ss in range(sss+1):
                            if h > -0.5 and (nn >= 1 or h > 0.5):
                                n = [c, 0, 2*h+nn, o, 0, 0, nn, 0, 1, 0, Cll, 0, 0, ss, 0, 0, 0, 0, 0]
                                mass = sum([n[i] * x for i, x in enumerate([C12, C13, H1, O16, O17, O18, N14, N15, H_pos, e_neg, Cl35, Cl37, F, S32O4, S34O4, Si28, Si29, Si30, I])])
                                lib.append([mass] + n)

                                if cal:
                                    n[8], n[9] = 0, 1  # Change to charge transfer mode
                                    mass = sum([n[i] * x for i, x in enumerate([C12, C13, H1, O16, O17, O18, N14, N15, H_pos, e_neg, Cl35, Cl37, F, S32O4, S34O4, Si28, Si29, Si30, I])])
                                    lib.append([mass] + n)

                                if extended:
                                    n[0], n[1] = c-1, 1  # 13C isotope
                                    mass = sum([n[i] * x for i, x in enumerate([C12, C13, H1, O16, O17, O18, N14, N15, H_pos, e_neg, Cl35, Cl37, F, S32O4, S34O4, Si28, Si29, Si30, I])])
                                    lib.append([mass] + n)

    # Sort the library by mass
    lib.sort(key=lambda x: x[0])

    return np.array(lib)
