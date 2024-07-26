import numpy as np

def averaging(data, masslist, engdata, engnames, ind, carryover=None):
    """
    Perform averaging on the input data.
    
    # Lines 8599-8600: Function definition and input parameters
    """
    
    # Lines 8604-8607: Extract indices and set up data structures
    indexM = ind['m']
    IS = np.column_stack((ind['s1'], ind['s2'], ind['s3']))
    ISM = [ind['s1M'], ind['s2M'], ind['s3M']]
    FinalDataNames = ['tSTRTav', 'tENDav', 'tSTRTsamp', 'tENDsamp', 'indM', 'I37/prim', 'I32/prim', 'I30/I32', 'prim_i/med(prim)', 'indS1', 'indS2', 'indS3', 'k-counter', 'k', '# avg'] + list(engnames) + [f"{m:.3f}" for m in masslist]
    
    # Lines 8608-8609: Initialize output arrays
    FinalData = np.zeros((1, len(FinalDataNames)))
    FinalDataErr = np.zeros((1, len(FinalDataNames)))

    # Lines 8611-8621: Initialize sampling parameters
    SAMPstart = np.zeros(3)
    SAMPend = np.zeros(3)
    SAMPvalue = np.zeros(3)
    LastSvalue = np.zeros(3)

    if carryover is not None:
        SAMPstart = carryover[0:3]
        SAMPend = carryover[3:6]
        SAMPvalue = carryover[6:9]
        LastSvalue = carryover[9:12]

    timeSTARTsamp = SAMPstart[0]
    timeENDsamp = SAMPend[0]

    # Lines 8623-8627: Set up counters and time variables
    counter = 0
    lleng = len(indexM)
    fil = np.max(np.where(np.abs(np.diff(indexM)) > 50)[0])
    if fil == -1:
        fil = lleng
    time = engdata[:, 0].astype(float)

    # Lines 8630-8644: Calculate Q/A parameters
    I21 = data[np.argmin(np.abs(masslist - 21.022))]
    I37 = data[np.argmin(np.abs(masslist - 37.028))]
    I38 = data[np.argmin(np.abs(masslist - 38.033))]
    I19 = I21 * 487.0 / np.max(corrtr(21, np.min(time)))
    I37[I37 < I38 * 645] = I38[I37 < I38 * 645] * 645
    I37 = I37 / np.max(corrtr(37, np.min(time)))
    Prim = I19 + I37
    Prim[Prim < 1e-7] = np.min(Prim) + 1e-7
    I34 = data[np.argmin(np.abs(masslist - 33.994))]
    I32 = data[np.argmin(np.abs(masslist - 31.989))]
    I32[I32 > 6000] = I34[I32 > 6000] * 242
    I32 = I32 / np.max(corrtr(32, np.min(time)))
    I30 = data[np.argmin(np.abs(masslist - 29.997))] / np.max(corrtr(30, np.min(time)))
    
    F37 = I37 / Prim
    F32 = I32 / Prim
    F30 = I30 / I32
    prStab = Prim / np.median(Prim)

    # Lines 8649-8653: Calculate time-related variables
    avTime = ind['T'] / (24 * 3600)  # in days
    tstart = time[0]
    tend = time[lleng - 1]
    step = (tend - tstart) / lleng
    stepSEC = step * (24 * 3600)

    # Lines 8656-8658: Initialize arrays for sample data
    Send10 = np.zeros((1000, 3))
    Send20 = np.zeros((1000, 3))
    Sstrt = np.zeros((1000, 3))

    # Lines 8660-8724: Main processing loop
    # ... (The rest of the function would be implemented here, following the IDL logic)

    # Return the results as a dictionary (equivalent to IDL structure)
    return {
        'names': FinalDataNames,
        'data': FinalData,
        'dataErr': FinalDataErr,
        'restdata': restdata,
        'resteng': resteng,
        'carryover': carryover
    }


def averaging_old(data, masslist, engdata, engnames, ind, carryover):
    # Line 8726-8729: Initialize variables
    indexM = ind.m
    IS = np.column_stack((ind.s1, ind.s2, ind.s3))
    ISM = [ind.s1M, ind.s2M, ind.s3M]
    FinalDataNames = ['tSTRTav', 'tENDav', 'tSTRTsamp', 'tENDsamp', 'indM', 'indS1', 'indS2', 'indS3', 'k-counter', 'k', '# avg'] + list(engnames) + [f"{m:.3f}" for m in masslist]
    
    # Line 8730-8731: Initialize FinalData and FinalDataErr arrays
    FinalData = np.zeros((1, len(FinalDataNames)))
    FinalDataErr = np.zeros((1, len(FinalDataNames)))

    # Line 8733-8736: Initialize sample arrays
    SAMPstart = np.zeros(3)
    SAMPend = np.zeros(3)
    SAMPvalue = np.zeros(3)
    LastSvalue = np.zeros(3)

    # Line 8738-8743: Handle carryover data if it exists
    if carryover is not None:
        SAMPstart = carryover[0:3]
        SAMPend = carryover[3:6]
        SAMPvalue = carryover[6:9]
        LastSvalue = carryover[9:12]

    # Line 8745-8746: Initialize sample start and end times
    timeSTARTsamp = SAMPstart[0]
    timeENDsamp = SAMPend[0]

    # Line 8748-8757: Calculate various time-related variables
    counter = 0
    lleng = len(indexM)
    fil = np.max(np.where(np.abs(np.diff(indexM)) > 50)[0]) + 1 if np.any(np.abs(np.diff(indexM)) > 50) else lleng
    time = engdata[:, 0].astype(float)
    avTime = float(ind.T) / (24.0 * 3600.0)  # in days
    tstart = time[0]
    tend = time[lleng - 1]
    step = (tend - tstart) / lleng
    stepSEC = step * (24.0 * 3600.0)

    # Line 8760-8762: Initialize arrays (not used in this function, can be removed)
    Send10 = np.zeros((1000, 3))
    Send20 = np.zeros((1000, 3))
    Sstrt = np.zeros((1000, 3))

    # Line 8764-8817: Main processing loop
    for k in range(min(fil + 2, lleng - 2)):
        if indexM[k] == indexM[k + 1] and (time[k] - tstart) < avTime:
            counter += 1
        else:
            if counter >= 1:
                timeSTARTav = time[k - counter]
                timeENDav = time[k]
                indM = np.mean(indexM[k - counter:k + 1])
                indS = np.zeros(3)
                
                for gg in range(3):
                    indxS = IS[:, gg]
                    if np.mean(indxS[k - counter:k + 1]) != LastSvalue[gg]:
                        filti = np.where(np.diff(np.append(indxS[k - counter:k + 1], indxS[k])) < 0)[0]
                        if indxS[k - counter] > LastSvalue[gg]:
                            SAMPstart[gg] = time[k - counter]
                        elif len(filti) > 0:
                            SAMPstart[gg] = time[filti[-1] + 1]
                        
                        filti = np.where(np.diff(np.append(indxS[k - counter:k + 1], indxS[k])) > 0)[0]
                        if indxS[k - counter] < LastSvalue[gg]:
                            SAMPend[gg] = time[k - counter]
                            SAMPvalue[gg] = LastSvalue[gg]
                        elif len(filti) > 0:
                            SAMPend[gg] = time[filti[-1]]
                            SAMPvalue[gg] = np.max(np.abs(indxS[k - counter:k + 1]))
                        
                        LastSvalue[gg] = indxS[k]

                    if 0 <= indM - ISM[gg] < 50:
                        timeSTARTsamp = SAMPstart[gg]
                        timeENDsamp = SAMPend[gg]
                        indS[gg] = SAMPvalue[gg]

                subset = np.concatenate(([timeSTARTav, timeENDav, timeSTARTsamp, timeENDsamp, indM], indS, 
                                         [k - counter, k, counter + 1], np.mean(engdata[k - counter:k + 1, :], axis=0)))
                FinalData = np.vstack((FinalData, np.concatenate((subset, np.mean(data[:, k - counter:k + 1], axis=1)))))
                FinalDataErr = np.vstack((FinalDataErr, np.concatenate((subset, np.sqrt(stepSEC * np.sum(data[:, k - counter:k + 1], axis=1)) / (stepSEC * np.sum(data[:, k - counter:k + 1], axis=1))))))
                
                last = k
            
            counter = 0
            tstart = time[k]

    # Line 8818-8826: Prepare return data
    carryover = np.concatenate((SAMPstart, SAMPend, SAMPvalue, LastSvalue))
    last = -1 if 'last' not in locals() else last
    restind = np.arange(last + 1, lleng)
    restdata = data[:, restind]
    resteng = engdata[restind, :]
    
    if FinalData.shape[0] > 1:
        FinalData = FinalData[1:, :]
        FinalDataErr = FinalDataErr[1:, :]

    return {
        'names': FinalDataNames,
        'data': FinalData,
        'dataErr': FinalDataErr,
        'restdata': restdata,
        'resteng': resteng,
        'carryover': carryover
    }

def cal3pt(peaklist, a, t0, ex, sampint, lib, mode, instrument):
    """
    3-point recalibration function
    """
    # Lines 8827-8834: Select appropriate library masses based on mode
    flib = -1
    if mode == 0:
        flib = np.where((2 * np.floor(np.floor(lib + 0.4) / 2) != np.floor(lib + 0.4)) | (lib < 50))[0]
        if flib.size > 0:
            lib = lib[flib]
    else:
        flib = np.where((2 * np.floor(np.floor(lib + 0.4) / 2) == np.floor(lib + 0.4)) | (lib < 50))[0]
        if flib.size > 0:
            lib = lib[flib]

    # Line 8839: Calculate massesFINE
    massesFINE = m2t(peaklist[:, 7], a, t0, ex, sampint, time=True)
    
    # Lines 8840-8842: Initialize variables
    a3, t03, ex3 = a, t0, ex

    # Lines 8844-8851: Get calibration masses from parameters
    M1a, M2a, M3a = getpar('M1a'), getpar('M2a'), getpar('M3a')
    M1b = M1a if getpar('ForceThroughMaSet') == 1 else getpar('M1b')
    M2b = M2a if getpar('ForceThroughMaSet') == 1 else getpar('M2b')
    M3b = M3a if getpar('ForceThroughMaSet') == 1 else getpar('M3b')
    Mset = np.array([M1a, M2a, M3a])
    tolli = getpar('tol_ppm')

    # Lines 8853-8855: Get ion parameters
    ions = np.array([getpar('ION1'), getpar('ION2'), getpar('ION3')])
    devs = np.array([getpar('DEV1'), getpar('DEV2'), getpar('DEV3')])

    # Lines 8858-8865: Initialize calibration matrix
    scores = 0
    Mcal = np.array([
        [M1a, M1a, M1a, M1a, M1b, M1b, M1b, M1b],
        [M2a, M2a, M2b, M2b, M2a, M2a, M2b, M2b],
        [M3a, M3b, M3a, M3b, M3a, M3b, M3a, M3b]
    ])

    # Lines 8867-8920: Main calibration loop
    for r in range(8):
        M1, M2, M3 = Mcal[:, r]
        
        I1 = np.where(np.abs(massesFINE - M1) < tolli * M1 * 1e-6)[0]
        I2 = np.where(np.abs(massesFINE - M2) < tolli * M2 * 1e-6)[0]
        I3 = np.where(np.abs(massesFINE - M3) < tolli * 3 * M3 * 1e-6)[0]

        if I1.size > 0 and I2.size > 0 and I3.size > 0:
            I1, I2, I3 = I1[0], I2[0], I3[0]
            T1 = peaklist[I1, 7] * sampint / 1e-10
            T2 = peaklist[I2, 7] * sampint / 1e-10
            T3 = peaklist[I3, 7] * sampint / 1e-10

            # Find optimal exponent
            Const = (T3 - T1) / (T2 - T1)
            ex4 = optimize_exponent(M1, M2, M3, Const)

            a4 = (T2 - T1) / (M2**ex4 - M1**ex4)
            t04 = T1 - a4 * M1**ex4

            if getpar('exMin') < ex4 < getpar('exMax'):
                tst = test_scale(m2t(peaklist[:, 7], a4, t04, ex4, sampint, time=True), lib, instrument, ions, devs)
                
                print([M1, M2, M3, tst['scorppm']])
                
                if tst['scorppm'] > scores:
                    scores, ex3, t03, a3, Mset = tst['scorppm'], ex4, t04, a4, [M1, M2, M3]

    # Lines 8922-8990: Additional optimization if not forced through MaSet
    if getpar('ForceThroughMaSet') != 1:
        # Variation of 3-point calibration
        varset_ppm = [20, 9, 3]
        if instrument.startswith('TOF1000'):
            varset_ppm = [x * 5 for x in varset_ppm]

        for s in range(3):
            var_ppm = varset_ppm[s]
            M1, M2, M3 = Mset
            
            # Generate variations for M1, M2, M3
            M1var, M2var, M3var = generate_mass_variations(M1, M2, M3, var_ppm)

            for r in range(64):
                M1, M2, M3 = M1var[r], M2var[r], M3var[r]
                
                # Perform calibration and testing similar to the main loop
                # ... (similar to the code in the main loop)

    return {'a': a3, 't0': t03, 'ex': ex3, 'Mset': Mset}


def calcppb(cpsdata, masslist, Prim1data, Prim2data, pdrift, udrift, tdrift, zeit=0):
    """
    Calculates volume mixing ratio
    """
    # Line 8997: Default zeit value
    if not np.any(zeit):
        zeit = 0

    # Lines 8999-9004: Reshape input data
    ppbdata = cpsdata.copy()
    Prim1 = np.reshape(Prim1data, -1)
    Prim2 = np.reshape(Prim2data, -1)
    pdrift = np.reshape(pdrift, -1)
    udrift = np.reshape(udrift, -1)
    tdrift = np.reshape(tdrift, -1)

    # Lines 9007-9008: Get parameters
    k19 = max(getpar('k19'))  # units e-9
    k37 = max(getpar('k37'))  # units e-9

    # Line 9010: Calculate transmission efficiency
    transe = corrtr(np.arange(2000), zeit)

    # Lines 9014-9015: Get dimensions of input data
    nDim = cpsdata.ndim
    AnzMass = len(masslist)

    # Lines 9017-9020: Handle 1D input cases
    if nDim == 1 and AnzMass == 1:
        ppbdata = np.reshape(ppbdata, -1)
        ppbdata = np.transpose([ppbdata, ppbdata])

    # Lines 9022-9025: Determine number of cycles
    if nDim == 1 and AnzMass > 1:
        cyc = 1
    elif np.max(ppbdata.shape) == 1:
        cyc = 1
    else:
        Dim = ppbdata.shape
        Cyc = Dim[1]

    # Lines 9027-9028: Get reaction parameters
    d = max(getpar('reactionlength'))  # reactionlength in cm
    mu0 = max(getpar('reduced_mobility'))  # reduced mobility of H3O+ in N2

    # Lines 9031-9033: Calculate reaction parameters
    mu = mu0 * (1013.25 / pdrift) * (tdrift / 273.15)  # [=] cm2/Vs
    trxn = d / ((mu * udrift) / d)  # reaction time in s
    Nmolec = 24.63 * 298 * pdrift / (1013.25 * tdrift)

    # Lines 9036-9040: Calculate ppb for each mass
    for i in range(AnzMass):
        aa = ppbdata[i, :] / transe[int(masslist[i])]
        bb = Prim1 * k19 / transe[21] + Prim2 * k37 / transe[37]
        cc = trxn * Nmolec
        ppbdata[i, :] = aa / (bb * cc)

    # Lines 9042-9043: Handle 1D output case
    if nDim == 1:
        ppbdata = ppbdata[0, :]

    # Line 9045: Return result
    return ppbdata

def calcppbold(cpsdata, masslist, Prim1data, Prim2data, pdrift, udrift, tdrift):
    """
    Calculates volume mixing ratio
    """
    # Line 9050-9055: Initialize variables
    ppbdata = cpsdata.copy()  # Create a copy to avoid modifying the original data
    Prim1 = Prim1data.ravel()  # Equivalent to IDL's reform()
    Prim2 = Prim2data.ravel()
    pdrift = pdrift.ravel()
    udrift = udrift.ravel()
    tdrift = tdrift.ravel()

    # Line 9058-9059: Get parameters
    k19 = max(getpar('k19'))  # units e-9
    k37 = max(getpar('k37'))  # units e-9

    # Line 9065-9066: Get dimensions of input data
    nDim = cpsdata.ndim
    AnzMass = max(masslist.shape)

    # Line 9068-9071: Handle 1D input cases
    if nDim == 1 and AnzMass == 1:
        ppbdata = ppbdata.reshape(-1)
        ppbdata = ppbdata[:, None].repeat(2, axis=1)

    # Line 9073-9076: Determine number of cycles
    if nDim == 1 and AnzMass > 1:
        cyc = 1
    elif max(ppbdata.shape) == 1:
        cyc = 1
    else:
        Dim = ppbdata.shape
        Cyc = Dim[1]

    # Line 9078-9079: Get reaction parameters
    d = max(getpar('reactionlength'))  # reactionlength in cm
    mu0 = max(getpar('reduced_mobility'))  # reduced mobility of H3O+ in N2

    # Line 9082-9084: Calculate reaction parameters
    mu = mu0 * (1013.25 / pdrift) * (tdrift / 273.15)  # [=] cm2/Vs
    trxn = d / ((mu * udrift) / d)  # reaction time in s
    Nmolec = 24.63 * 298 * pdrift / (1013.25 * tdrift)

    # Line 9087-9089: Calculate ppb values
    for i in range(AnzMass):
        ppbdata[i, :] = ((ppbdata[i, :] / corrtr(masslist[i])) / 
                         (Prim1 * k19 / corrtr(21) + Prim2 * k37 / corrtr(38))) / (trxn * Nmolec)

    # Line 9091: Handle 1D output case
    if nDim == 1:
        ppbdata = ppbdata[0, :]

    # Line 9094: Return result
    return ppbdata

def cal_crude(peaklist, sampint, lib, instrument):
    """
    Crude/Rough Mass scale calibration
    """
    # Initialize variables (lines 9103-9105)
    a = -9999
    t0 = -9999
    mean_dev = 1000
    
    # Get number of peaks to consider (line 9106)
    nnn = max(getpar('PeaksToConsider'))
    
    # Sort peaklist and get top peaks (lines 9107-9110)
    dims = peaklist.shape
    top16 = peaklist[peaklist[:, 6].argsort()]
    top16all = top16[dims[0]-nnn:dims[0], :]
    top16 = top16all[:, 7]
    
    # Set tolerance values (lines 9111-9116)
    tol_crude_low = 0.02
    tol_crude_high = 0.1
    if instrument.startswith('TOF1000'):
        tol_crude_low = 0.1
        tol_crude_high = 0.3
    
    # Define mass pairs (lines 9117-9118)
    crude_pairs = [[19.018, 37.028], [29.997, 45.992], [19.018, 29.997], 
                   [31.989, 45.992], [19.018, 31.989], [37.028, 59.049]]
    cal_pairs = [[21.022, 59.049], [30.994, 47.997], [21.022, 30.994], 
                 [33.994, 47.997], [21.022, 33.994], [39.033, 60.052]]
    
    print('MassLow, timeLow, MassHigh, timeHigh, mean deviation, a, t0, adopted')
    
    # Main calibration loop (lines 9124-9170)
    ii = 0
    for i in range(nnn):
        time_low = max(float(top16[i]))
        for j in range(nnn):
            time_high = max(float(top16[j]))
            
            if time_high > time_low:
                for k in range(len(crude_pairs)):
                    mass_low, mass_high = crude_pairs[k]
                    mass_low_cal, mass_high_cal = cal_pairs[k]
                    
                    # Calculate t0 and a (lines 9136-9137)
                    tt0 = (mass_low**0.5 * time_high * (sampint/1e-10) - mass_high**0.5 * time_low * (sampint/1e-10)) / (mass_low**0.5 - mass_high**0.5)
                    aa = (time_low * (sampint/1e-10) - tt0) / mass_low**0.5
                    
                    # Convert peaklist times to masses (line 9138)
                    masses = m2t(peaklist[:, 7], aa, tt0, 0.5, sampint, time=True)
                    
                    # Check if masses are within tolerance (line 9141)
                    if min(abs(masses - mass_low_cal)) < tol_crude_low and min(abs(masses - mass_high_cal)) < tol_crude_high:
                        # Refine calibration (lines 9144-9151)
                        i1 = abs(masses - mass_low_cal).argsort()
                        time_low_cal = peaklist[i1[0], 7]
                        i2 = abs(masses - mass_high_cal).argsort()
                        time_high_cal = peaklist[i2[0], 7]
                        
                        tt0 = max((mass_low_cal**0.5 * time_high_cal * (sampint/1e-10) - mass_high_cal**0.5 * time_low_cal * (sampint/1e-10)) / (mass_low_cal**0.5 - mass_high_cal**0.5))
                        aa = max((time_low_cal * (sampint/1e-10) - tt0) / mass_low_cal**0.5)
                        
                        masses = m2t(peaklist[:, 7], aa, tt0, 0.5, sampint, time=True)
                        test = abs(testscale3(masses, lib, instrument))
                        print(test, mass_low_cal, mass_high_cal, aa, tt0)
                        
                        # Update best calibration if test is better (lines 9155-9161)
                        if test < mean_dev:
                            mean_dev = test
                            a = aa
                            t0 = tt0
                        
                        ii += 1
    
    # Print final results (lines 9175-9178)
    print('____________________________')
    print('mean deviation, a, t0')
    print([mean_dev, a, t0])
    print('____________________________')
    
    # Return results as a dictionary (line 9182)
    return {'a': a, 't0': t0, 'ex': 0.5, 'Top16': top16all}

def cal_fine(peaktimes, lib, a, t0, samp_int, mode, instrument):
    # Line 9187: Initialize variables
    analytics = 0
    rat120 = -9999
    ex = 0.5

    # Lines 9191-9192: Get ion parameters
    ions = [getpar('ION1'), getpar('ION2'), getpar('ION3')]
    devs = [getpar('DEV1'), getpar('DEV2'), getpar('DEV3')]

    # Line 9196: Outer loop for refinement
    for r in range(4):
        a_raw, t0_raw, ex_raw = a, t0, ex
        
        # Line 9201: Inner loop for parameter adjustment
        for q in range(-3, 4):
            mag = [1, 3.5, 12, 43]
            ex1 = ex_raw + 0.00001 * q / mag[r]
            grid = 29
            stepa = 50.002 / (grid * mag[r])
            stept0 = 500.01 / (grid * mag[r])
            
            # Line 9208: Ensure grid is odd
            if grid % 2 == 0:
                grid += 1
            
            # Lines 9209-9210: Initialize score arrays
            sco = [[0 for _ in range(grid)] for _ in range(grid)]
            sco120 = [[0 for _ in range(grid)] for _ in range(grid)]

            # Line 9212: Loop over grid
            for k in range(-grid//2, grid//2 + 1):
                if ex1 > max(getpar('exMin')) and ex1 < max(getpar('exMax')):
                    aa = a_raw + stepa * k
                    
                    # Lines 9217-9222: Calculate scores
                    for l in range(-grid//2, grid//2 + 1):
                        t00 = t0_raw + stept0 * l
                        tst = test_scale2(m2t(peaktimes, aa, t00, ex1, samp_int, time=True), lib, instrument, ions, devs)
                        sco[k + grid//2][l + grid//2] = tst.scor
                        sco120[k + grid//2][l + grid//2] = tst.scor120

            # Lines 9226-9234: Find maximum score
            scomist = [row[:] for row in sco]
            maxi = max(range(len(scomist)), key=lambda i: max(scomist[i]))
            maxx120 = max(max(row) for row in sco120)
            
            while sco120[maxi // grid][maxi % grid] < 0.9 * maxx120:
                scomist[maxi // grid][maxi % grid] = 0
                maxi = max(range(len(scomist)), key=lambda i: max(scomist[i]))

            # Lines 9235-9241: Calculate new parameters
            l, k = divmod(maxi, grid)
            a1 = max(a_raw + stepa * (k - grid//2))
            t01 = max(t0_raw + stept0 * (l - grid//2))
            old = test_scale2(m2t(peaktimes, a, t0, ex, samp_int, time=True), lib, instrument, ions, devs)
            new = test_scale2(m2t(peaktimes, a1, t01, ex1, samp_int, time=True), lib, instrument, ions, devs)

            # Lines 9244-9250: Update parameters if score improved
            if new.scor > old.scor and new.scor120 > old.scor120 * 0.9:
                a, t0, ex = a1, t01, ex1
                rat120 = new.rat120

            # Lines 9254-9314: Analytics section (omitted for brevity)

    # Lines 9318-9319: Return results
    return {'a': a, 't0': t0, 'ex': ex, 'rat120': rat120}


def corr_poiss_dead(tof_data, extractions, sampint):
    # Line 9325-9326: Get dimensions of input data
    dim = tof_data.shape
    bins, cyc = dim[0], dim[1]
    
    # Line 9327-9328: Calculate non-extending and extending dead times
    nedt = int(getpar('NonExtendingDeadTime') * 10 * 1e-10 / sampint)
    edt = int(getpar('ExtendingDeadTime') * 10 * 1e-10 / sampint)
    
    # Line 9329: Loop over cycles
    for rr in range(1, cyc):
        tof_data_cor = tof_data[:, rr].copy()
        
        # Line 9331: Find indices where count rate exceeds threshold
        indx = np.where(1000 * 1e-10 * tof_data[:, rr] / (extractions * sampint) > 1)[0]
        
        # Line 9332-9333: Loop over selected indices
        for ss in range(len(indx)):
            si = tof_data[indx[ss], rr]
            r = extractions
            
            # Line 9336-9337: Calculate Sj and nj
            sj = np.sum(tof_data[indx[ss]-nedt:indx[ss]-edt, rr]) / r
            nj = np.sum(tof_data_cor[indx[ss]-edt:indx[ss]]) / r
            
            # Line 9338: Calculate corrected count (ni)
            ni = -r * np.log(1 - (np.exp(nj) * si / (r * (1 - sj))))
            
            # Line 9344: Update corrected data
            tof_data_cor[indx[ss]] = ni
        
        # Line 9346: Update original data with corrected values
        tof_data[:, rr] = tof_data_cor
    
    # Line 9349: Return corrected data
    return tof_data

def corrtr(mass, zeit):
    """
    Calculate the mass transmission of the TOF from m0-m1500.
    
    Args:
    mass (float or array): Mass values
    zeit (float): Time value
    
    Returns:
    array: Transmission values
    """
    # Line 9354: Convert mass to float
    mass = float(mass)

    # Line 9357: Get transmission parameters for the current instrument and time
    pari = gettr(getInstID(), zeit)

    # Lines 9358-9364: Handle 5-parameter case
    if len(pari) == 5:
        ex, mLOW, wLOW, mHigh, wHigh = pari
        tran = tr(mass, ex, mLOW, wLOW, mHigh, wHigh, norm=True)
    
    # Lines 9365-9374: Handle 13-parameter case
    elif len(pari) == 13:
        tranL = tranH = -1
        low = mass <= pari[7]
        high = mass > pari[7]

        if any(low):
            tranL = (pari[6] + pari[5]*mass[low] + pari[4]*mass[low]**2 + 
                     pari[3]*mass[low]**3 + pari[2]*mass[low]**4 + 
                     pari[1]*mass[low]**5 + pari[0]*mass[low]**6)
        
        if any(high):
            tranH = (pari[12] + pari[11]*mass[high] + pari[10]*mass[high]**2 + 
                     pari[9]*mass[high]**3 + pari[8]*mass[high]**4)
        
        tran = [tranL, tranH]
        tran = [t for t in tran if t != -1]

    # Line 9379: Return the calculated transmission values
    return tran

def corrtr_untilDEC2020(mass, zeit):
    # Line 9384: Convert mass to float
    mass = float(mass)
    
    # Lines 9385-9387: Check if UniqueCampaignName is default
    if 'default' in getstr('UniqueCampaignName'):
        # Calculate transmission using polynomial coefficients
        tran = (max(getpar('P0')) + max(getpar('P1'))*mass + max(getpar('P2'))*mass**2 +
                max(getpar('P3'))*mass**3 + max(getpar('P4'))*mass**4 + max(getpar('P5'))*mass**5 +
                max(getpar('P6'))*mass**6)
    else:
        # Lines 9389-9390: Get transmission parameters for the campaign
        pari = gettr(getstr('UniqueCampaignName'), zeit)
        
        # Lines 9391-9397: Handle case with 5 parameters
        if len(pari) == 5:
            ex, mLOW, wLOW, mHigh, wHigh = pari
            tran = tr(mass, ex, mLOW, wLOW, mHigh, wHigh, norm=True)
        
        # Lines 9398-9408: Handle case with 13 parameters
        elif len(pari) == 13:
            tranL = tranH = -1
            low = mass <= pari[7]
            high = mass > pari[7]
            
            if any(low):
                tranL = (pari[6] + pari[5]*mass[low] + pari[4]*mass[low]**2 + pari[3]*mass[low]**3 +
                         pari[2]*mass[low]**4 + pari[1]*mass[low]**5 + pari[0]*mass[low]**6)
            
            if any(high):
                tranH = (pari[12] + pari[11]*mass[high] + pari[10]*mass[high]**2 +
                         pari[9]*mass[high]**3 + pari[8]*mass[high]**4)
            
            tran = [tranL, tranH]
            tran = [t for t in tran if t != -1]
    
    # Line 9411: Return the calculated transmission
    return tran

def corrtr_new(mass):
    # Line 9414: Convert mass to float
    mass = float(mass)
    
    # Line 9417: Call the tr function with specific parameters
    # Note: The '/norm' keyword in IDL is likely handled within the tr function in Python
    tr_value = tr(mass, 0.6, 0, 0, 599.5, 299)
    
    # Line 9418: Return the calculated transmission value
    return tr_value

# Lines 9413 and 9419: Function definition and end
# The 'end' statement in IDL is not needed in Python as function blocks are defined by indentation


def tr(mass, ex, MASSlow, WIDTHlow, MASShigh, WIDTHhigh, norm=False):
    """
    Calculate the mass transmission of the TOF from m0-m1500.
    
    Line 9423: Function definition with parameters matching IDL version
    """
    # Line 9424-9428: Convert inputs to float
    mass = float(mass)
    MASSlow = float(MASSlow)
    WIDTHlow = float(WIDTHlow)
    MASShigh = float(MASShigh)
    WIDTHhigh = float(WIDTHhigh)

    # Line 9429: Determine offset based on WIDTHlow
    offset = 0 if WIDTHlow >= 0 else 1

    # Line 9430: Calculate transmission
    tr_value = (mass**ex * 
                (offset + 1.0 / (1.0 + exp(-(mass - MASSlow) / WIDTHlow))) * 
                (1.0 / (1.0 + exp((mass - MASShigh) / WIDTHhigh))))

    # Line 9431-9432: Normalize if norm is True
    if norm:
        tr_value /= (59**ex * 
                     (offset + 1.0 / (1.0 + exp(-(59 - MASSlow) / WIDTHlow))) * 
                     (1.0 / (1.0 + exp((59 - MASShigh) / WIDTHhigh))))

    # Line 9434: Return the calculated value
    return tr_value

def detect_peaks(sum_spectrum, a, t0, ex, samp_int, instrument, clean=False):
    # Line 9437-9438: Function definition and optional 'clean' parameter
    
    # Line 9441-9445: Smooth the sum spectrum and get derivative data
    smoo = smooth_sum_spec(sum_spectrum, a, t0, ex, samp_int, instrument)
    data2 = smoo['data']
    deriv_data = smoo['deriv_data']
    data3 = deriv(deriv_data)
    n_pts = smoo['n_pts']
    
    # Line 9448-9453: Initialize variables for peak detection
    start = max(floor(m2t(get_par('Min_Mass'), a, t0, ex, samp_int)), 1000)
    ende = n_pts - 1
    peak_list = [[0] * 9 for _ in range(20000)]
    peak_count = 0
    i = 1 + start + 200
    run_num = 0
    
    # Line 9454-9552: Main peak detection loop
    while i < ende - 1000:
        # Calculate mass unit and step size
        m_unit = floor(i - m2t(m2t(i, a, t0, ex, samp_int, time=True) - 1, a, t0, ex, samp_int, mass=True))
        step = max(1, floor(m_unit / 15))
        
        # Calculate thresholds for peak detection
        per1 = deriv_data[max(0, i - m_unit):i]
        f1 = sum_spectrum[max(0, i - m_unit):i]
        per2 = deriv_data[i:min(ende, i + m_unit)]
        f2 = sum_spectrum[i:min(ende, i + m_unit)]
        thresh_a = median(abs(per1[f1 != 0])) if any(f1 != 0) else 0
        thresh_b = median(abs(per2[f2 != 0])) if any(f2 != 0) else 0
        thresh = min(thresh_a, thresh_b) * 6
        
        i += step
        condi = max(deriv_data[i:i+step]) > thresh
        
        if condi:
            run_num += 1
            # Find peak start, max slope, peak max, and peak end
            i_peak_start, i_slope_max, i_peak_max, i_slope_min, i_peak_end = find_peak_characteristics(deriv_data, i, step, thresh)
            
            # Calculate peak properties
            peak_broadness = m2t(i_peak_max, a, t0, ex, samp_int, time=True) / (m2t(i_peak_end, a, t0, ex, samp_int, time=True) - m2t(i_peak_start, a, t0, ex, samp_int, time=True))
            slope_max = deriv_data[i_slope_max]
            slope_min = deriv_data[i_slope_min]
            counts_max = data2[i_peak_max]
            der2 = data3[i_peak_max]
            counts_min = min(data2[i_peak_start:i_peak_end])
            
            # Apply cleaning if required
            if clean:
                good = clean_peak(data2, data3, i_peak_max, m_unit, peak_broadness, counts_max)
            else:
                good = True
            
            # Store peak information if it's a good peak
            if good:
                peak_list[peak_count] = [
                    m2t(i_peak_max, a, t0, ex, samp_int, time=True),
                    i_peak_start, i_peak_end, peak_broadness,
                    slope_max, slope_min, counts_max, i_peak_max, counts_min
                ]
                peak_count += 1
                i = i_peak_end - step + 1
            elif i_peak_max > 1:
                i = i_peak_max - step + 1
    
    # Line 9553: Return the list of detected peaks
    return peak_list[:peak_count]

def index(ppbdata, masses2, engnms, entdt, indexfile):
    # Line 9557: Read index file
    IndexPAR = ReadIndexFile(IndexFile)
    
    # Line 9560: Trim whitespace from engineering names
    engnms = [name.strip() for name in engnms]
    
    # Lines 9561-9567: Initialize index arrays
    INDM = [0] * len(entdt)
    INDS1 = [0] * len(entdt)
    INDS2 = [0] * len(entdt)
    INDS3 = [0] * len(entdt)
    INDS1M = 0
    INDS2M = 0
    INDS3M = 0
    
    # Line 9570: Get vector index from IndexPAR
    Dind = IndexPAR.VECTORIND
    
    # Lines 9571-9576: Handle default case
    if max(str2vec(Dind[0])) == -9999:
        INDM = [300] * len(INDM)
        MaxAvTime = STR2VEC(IndexPAR.INTERVAL) if 'INTERVAL' in IndexPAR.__dict__ else -9999
    else:
        # Lines 9577-9581: Get parameters from IndexPAR
        OPCO = IndexPAR.CONDIND
        VAL = IndexPAR.VALUEIND
        IVAL = IndexPAR.IFIND
        IDEF = IndexPAR.ELSEIND
        MaxAvTime = STR2VEC(IndexPAR.INTERVAL)
        
        # Lines 9583-9628: Construct INDM
        for i in range(10):
            if 0.5 < max(str2vec(OPCO[i])) < 5.5:
                indi = str2vec(Dind[i])
                n = len(indi)
                
                # Handle data selection based on indi values
                if max(indi[0]) * 10000 - int(max(indi[0])) * 10000 > 0:
                    data = sum(ppbdata[abs(masses2 - max(indi[j])) < 0.001] for j in range(n))
                else:
                    data = sum(entdt[:, indi[j]] for j in range(n))
                
                # Special handling for time data
                if max(indi) == 0 and n == 1:
                    sttrt = getpar('StartExperiment')
                    perr = getpar('PeriodExperiment')
                    data = (data - sttrt) * 24 * 60
                    data = data % perr
                
                # Apply filters based on OPCO values
                if str2vec(idef[i]) != -9999:
                    INDM = [str2vec(idef[i])] * len(INDM)
                
                filter_conditions = {
                    1: lambda x: x < max(str2vec(VAL[i])),
                    2: lambda x: x == max(str2vec(VAL[i])),
                    3: lambda x: x > max(str2vec(VAL[i])),
                    4: lambda x: x == n * max(str2vec(VAL[i])),
                    5: lambda x: x >= max(str2vec(VAL[i]))
                }
                
                filter_condition = filter_conditions.get(max(str2vec(OPCO[i])))
                if filter_condition:
                    filter = [idx for idx, val in enumerate(data) if filter_condition(val)]
                
                # Apply additional filtering for specific cases
                if len(filter) > 9 and max(indi[0]) * 10000 - int(max(indi[0])) * 10000 > 0:
                    tol = 9
                    if len(filter) > tol and max(indi[0]) * 10000 - int(max(indi[0])) * 10000 > 0:
                        if filter[0 + tol] > filter[0] + tol + 1:
                            filter[0] = -1
                        for zz in range(len(filter)):
                            if filter[min(zz + tol, len(filter) - 1)] > filter[zz] + tol + 1:
                                if filter[max(0, zz - tol)] < filter[zz] - tol - 1:
                                    filter[zz] = -1
                        if filter[len(filter) - 1 - tol] < filter[len(filter) - 1] - tol - 1:
                            filter[len(filter) - 1] = -1
                        filter = [f for f in filter if f > -0.5]
                        filter.extend([f - 1 for f in filter if filter[filter.index(f) + 1] - f == 1])
                        filter.extend([f - 1 for f in filter if filter[filter.index(f) + 1] - f == 2])
                        filter.extend([f - 2 for f in filter if filter[filter.index(f) + 1] - f == 2])

                   
                
                if filter:
                    for idx in filter:
                        INDM[idx] = str2vec(IVAL[i])
            
            # Lines 9629-9639: Handle additional conditions
            elif 5.5 < max(str2vec(OPCO[i])) < 8.5:
                indi = str2vec(Dind[i])
                data = entdt[:, indi[0]]
                n = len(indi)
                if n > 2:
                    for j in range(1, n-2):
                        data += entdt[:, indi[j]]
                if max(str2vec(OPCO[i])) < 7.5 and n >= 2:
                    data += abs(entdt[:, indi[n-1]] - 1)
                if max(str2vec(OPCO[i])) == 8 and n >= 2:
                    data -= entdt[:, indi[n-1]]
                if str2vec(idef[i]) != -9999:
                    INDM = [str2vec(idef[i])] * len(INDM)
                if max(str2vec(OPCO[i])) == 6:
                    filter = [idx for idx, val in enumerate(data) if val == n * max(str2vec(VAL[i]))]
                elif max(str2vec(OPCO[i])) == 7:
                    filter = [idx for idx, val in enumerate(data) if val >= max(str2vec(VAL[i]))]
                elif max(str2vec(OPCO[i])) == 8:
                    filter = [idx for idx, val in enumerate(data) if val > max(str2vec(VAL[i]))]
                if filter:
                    for idx in filter:
                        INDM[idx] = str2vec(IVAL[i])
        
        # Lines 9642-9678: Apply hardcopy to INDM
        INDVALUE = IndexPAR.INDVALUE
        VALUES2ADD = IndexPAR.VALUES2ADD
        length = len(INDM)

        for i in range(4):
            INDVAL2 = max(str2vec(INDVALUE[i]))
            if INDVAL2 != -9999:
                VALUES2ADD2 = str2vec(VALUES2ADD[i])
                if max(VALUES2ADD2[0]) == -1:
                    VALUES2ADD2 = VALUES2ADD2[1:]
                    VALUES2ADD2 = VALUES2ADD2 * 8  # Repeat 8 times

                if len(VALUES2ADD2) == 1:
                    add = entdt[:, indi[0]]
                    INDM = [val + add[idx] if val == INDVAL2 else val for idx, val in enumerate(INDM)]
                    filli = [idx for idx in range(1, length) if abs(INDM[idx] - INDM[idx-1]) > 0]
                    filli = [idx for idx in filli if abs(INDM[idx] - INDVAL2) < 50]
                    for idx in filli:
                        INDM[idx] = INDVAL2
                else:
                    indstart = [idx+1 for idx in range(length-1) if INDM[idx+1] == INDVAL2 and INDM[idx] != INDM[idx+1]]
                    if indstart or INDM[0] == INDVAL2:
                        if not indstart:
                            indstart = [1]
                        indend = [idx for idx in range(length-1) if INDM[idx] == INDVAL2 and INDM[idx] != INDM[idx+1]]
                        if indend:
                            if indstart[0] > indend[0]:
                                indstart = [0] + indstart
                            if indstart[-1] > indend[-1]:
                                indend.append(length - 1)
                        elif indstart:
                            indend = [length - 1]

                        add = VALUES2ADD2 + [0] * 900  # Extend with zeros
                        for start, end in zip(indstart, indend):
                            segment_length = end - start + 1
                            INDM[start:end+1] = [val + add[idx] for idx, val in enumerate(INDM[start:end+1])]

       # Lines 9680-9734: Construct INDS
        for i in range(10, 16):
            if 0.5 < max(str2vec(OPCO[i])) < 5.5:
                indi = str2vec(Dind[i])
                data = entdt[:, indi[0]]
                n = len(indi)
                if n > 1:
                    for j in range(1, n):
                        data += entdt[:, indi[j]]
                
                if i == 10 and str2vec(idef[i]) != -9999:
                    indS1M = str2vec(idef[i])
                elif i == 12 and str2vec(idef[i]) != -9999:
                    indS2M = str2vec(idef[i])
                elif i == 14 and str2vec(idef[i]) != -9999:
                    indS3M = str2vec(idef[i])
                
                filter_conditions = {
                    1: lambda x: x < max(str2vec(VAL[i])),
                    2: lambda x: x == max(str2vec(VAL[i])),
                    3: lambda x: x > max(str2vec(VAL[i])),
                    4: lambda x: x == n * max(str2vec(VAL[i])),
                    5: lambda x: x >= max(str2vec(VAL[i]))
                }
                
                filter_condition = filter_conditions.get(max(str2vec(OPCO[i])))
                if filter_condition:
                    filter = [idx for idx, val in enumerate(data) if filter_condition(val)]
                
                if i in [10, 11] and max(filter) > -0.5:
                    indS1[filter] = str2vec(IVAL[i])
                elif i in [12, 13] and max(filter) > -0.5:
                    indS2[filter] = str2vec(IVAL[i])
                elif i in [14, 15] and max(filter) > -0.5:
                    indS3[filter] = str2vec(IVAL[i])
            
            elif 5.5 < max(str2vec(OPCO[i])) < 8.5:
                indi = str2vec(Dind[i])
                data = entdt[:, indi[0]]
                n = len(indi)
                if n > 2:
                    for j in range(1, n-2):
                        data += entdt[:, indi[j]]
                
                if max(str2vec(OPCO[i])) < 7.5 and n >= 2:
                    data += abs(entdt[:, indi[n-1]] - 1)
                elif max(str2vec(OPCO[i])) == 8 and n >= 2:
                    data -= entdt[:, indi[n-1]]
                
                if i == 10 and str2vec(idef[i]) != -9999:
                    indS1M = str2vec(idef[i])
                elif i == 12 and str2vec(idef[i]) != -9999:
                    indS2M = str2vec(idef[i])
                elif i == 14 and str2vec(idef[i]) != -9999:
                    indS3M = str2vec(idef[i])
                
                filter_conditions = {
                    6: lambda x: x == n * max(str2vec(VAL[i])),
                    7: lambda x: x >= max(str2vec(VAL[i])),
                    8: lambda x: x > max(str2vec(VAL[i]))
                }
                
                filter_condition = filter_conditions.get(max(str2vec(OPCO[i])))
                if filter_condition:
                    filter = [idx for idx, val in enumerate(data) if filter_condition(val)]
                
                if i in [10, 11] and max(filter) > -0.5:
                    indS1[filter] = str2vec(IVAL[i])
                elif i in [12, 13] and max(filter) > -0.5:
                    indS2[filter] = str2vec(IVAL[i])
                elif i in [14, 15] and max(filter) > -0.5:
                    indS3[filter] = str2vec(IVAL[i])
 
    # Line 9735: Create and return the result structure
    s1 = {
        'm': INDM,
        's1': INDS1,
        's2': INDS2,
        's3': INDS3,
        's1M': INDS1M,
        's2M': INDS2M,
        's3M': INDS3M,
        'T': MaxAvTime
    }
    return s1

def DeDrift(data, a, t0, SampInt, rest, maxis, t1start, t2start, m1, m2, n_junks, cur_junk):
    # Line 9738: Set developing flag
    developing = 0

    # Lines 9740-9741: Set smoothing parameters
    smoothy = 9
    smoothyhalf = round(float(smoothy - 1) / 2)

    # Lines 9742-9743: Add rest data if provided
    if max(rest) != -9999:
        data = np.column_stack((rest, data))

    # Lines 9746-9748: Calculate dimensions and splits
    dims = data.shape
    summi = np.sum(data, axis=1)
    split = [0] + [round(m2t(m, a, t0, 0.5, SampInt)) for m in [40, 100, 250]] + [dims[0]]

    # Lines 9750-9752: Initialize arrays
    maxarr = np.zeros((dims[1], 100))
    maxarrval = np.zeros((dims[1], 100))
    array = np.zeros((13, 100))

    # Lines 9753-9854: Main processing loop
    for spl in range(4):
        indx = range(split[spl], split[spl+1])
        if max([i in indx for i in [round(t1start), round(t2start)]]) > -0.5 or t1start == -9999:
            dat2 = data[indx, :]
            dat3 = data[indx, smoothyhalf]
            for j in range(10):  # track the timebin-shift of the 200 largest peaks
                maxbin = np.argmax(dat3)
                dat3[max(maxbin-50):maxbin+50] = 0  # set +/- 50 timbins around highest signal to zero 0
                
                if t1start != -9999:
                    while min(abs(maxbin+split[spl]-t1start), abs(maxbin+split[spl]-t2start)) > 25:
                        maxbin = np.argmax(dat3)
                        dat3[max(maxbin-50):maxbin+50] = 0
                    j = 9
                
                if min(abs(maxbin+split[spl]-t1start), abs(maxbin+split[spl]-t2start)) < 25 or t1start == -9999:
                    cyc_dat = dat2[:, 0]
                    segment = cyc_dat[max(maxbin-50):maxbin+50]
                    dat2[max(maxbin-50):maxbin+50, 0] = 0
                    
                    maxval = max(cyc_dat[maxbin])  # value of highest signal in cycle 0
                    ttest = tm_test(segment[[47,48,49,51,52,53]], segment[[0,1,2,97,98,99]])
                    r, coeff, fff, ggg, hhh = gaussfit(range(-50, 50), segment, nterms=4)
                    xival = fff / maxval
                    yerrval = ggg / maxval
                    A1val = hhh[1] / coeff[1]
                    A2val = hhh[2] / coeff[2]
                    pval = ttest[1]
                    ppval = ttest[0]
                    maxbin += split[spl]
                    
                    if developing == 1:
                        # Plotting code here
                        pass
                    
                    for i in range(1, dims[1]):
                        cyc_dat = dat2[:, i]
                        segment = cyc_dat[max(maxbin[i-1]-split[spl]-50):min(maxbin[i-1]-split[spl]+50, len(cyc_dat))]
                        if developing == 1:
                            # More plotting code here
                            pass
                        addi = np.argmax(segment)
                        addi = addi[abs(addi-50) == min(abs(addi-50))]
                        if max(segment) <= 0:
                            addi = 50
                        maxi = maxbin[i-1] - split[spl] - 50 + addi[0]
                        dat2[max(maxi-50):min(maxi+50, dat2.shape[0]), i] = 0
                        vali = max(cyc_dat[maxi])
                        ttest = tm_test(segment[[47,48,49,51,52,53]], segment[[0,1,2,97,98,99]])
                        pval = np.append(pval, ttest[1])
                        ppval = np.append(ppval, ttest[0])
                        r, coeff, fff, ggg, hhh = gaussfit(range(-50, 50), segment, nterms=4)
                        xival = np.append(xival, fff/vali)
                        yerrval = np.append(yerrval, ggg/vali)
                        A1val = np.append(A1val, hhh[1]/coeff[1])
                        A2val = np.append(A2val, hhh[2]/coeff[2])
                        maxbin = np.append(maxbin, maxi+split[spl])
                        maxval = np.append(maxval, vali)
                    
                    if developing == 1:
                        # More plotting code here
                        pass
                    
                    array[:, j+spl*10] = [j+spl*10, len(maxbin), np.mean(maxbin), 1000.0*np.std(maxbin)/np.mean(maxbin),
                                        np.mean(maxval), max(A1val), max(A2val), 100*max(pval), min(ppval),
                                        100.0*np.percentile(yerrval, 90), np.percentile(xival, [80, 90]), np.median(maxval)]
                    maxarr[:, j+spl*10] = maxbin
                    maxarrval[:, j+spl*10] = maxval

    # Lines 9856-9858: Post-processing of array
    array = array[:, :40]
    array[array == 0] = 9999
    if developing == 1:
        print(array)

    # Lines 9860-9869: Filter good peaks
    filt = np.where((array[7, :] < 5) & (array[11, :] < 10))[0]
    if cur_junk != 1:
        filt = np.where(array[7, :] < 9998)[0]
        # Print warnings for potentially lost peaks

    # Lines 9872-9874: Apply filter to arrays
    array = array[:, filt]
    maxarr = maxarr[:, filt]
    maxarrval = maxarrval[:, filt]

    # Lines 9876-9885: Plotting (if developing)
    if developing == 1:
        import matplotlib.pyplot as plt

        x = (maxarr[:, 0] - np.median(maxarr[:, 0])) * 100 + np.median(maxarr[:, 0])
        plt.figure(figsize=(16, 5))
        plt.plot(x, range(900))
        plt.plot(smooth(x, 9), range(900), color='red', linewidth=3)
        
        for j in range(1, len(filt)):
            x = (maxarr[:, j] - np.median(maxarr[:, j])) * 100 + np.median(maxarr[:, j])
            plt.plot(x, range(900))
            plt.plot(smooth(x, 9), range(900), color='red', linewidth=3)
        
        plt.xlim(0, 50000)
        plt.title('Maxarr Plot')
        plt.xlabel('Adjusted Time Bin')
        plt.ylabel('Index')
        plt.show()

    # Lines 9892-9895: Calculate mean time bins and vectors
    meanTbin = np.mean(maxarr, axis=0)
    T1vec = maxarr[:, np.argmin(meanTbin)]
    T2vec = maxarr[:, np.argmax(meanTbin)]

    # Lines 9897-9904: More plotting (if developing)
    if developing == 1:
        x = (T1vec - np.median(T1vec)) * 100 + np.median(T1vec)
        plt.plot(x, range(900))
        plt.plot(smooth(x, 9), range(900), color='blue', linewidth=3)

        x = (T2vec - np.median(T2vec)) * 100 + np.median(T2vec)
        plt.plot(x, range(900))
        plt.plot(smooth(x, 9), range(900), color='green', linewidth=3)

        plt.title('T1vec and T2vec Plot')
        plt.xlabel('Adjusted Time Bin')
        plt.ylabel('Index')
        plt.show()

    # Lines 9910-9913: Set rest and start values
    i = len(T1vec) - smoothyhalf - 1
    rest = data[:, i-smoothyhalf:i+smoothyhalf+1]
    t1start = T1vec[i]
    t2start = T2vec[i]

    # Lines 9915-9919: Plotting (if developing)

    if developing == 1:
        plt.figure(figsize=(12, 8))
        plt.subplot(2, 1, 1)
        plt.plot(t1vec, range(len(t2vec)))
        plt.plot(smooth(t1vec, smoothy), range(len(t2vec)), color='red', linewidth=3)
        plt.title('T1vec Plot')
        plt.ylabel('Index')

        plt.subplot(2, 1, 2)
        plt.plot(t2vec, range(len(t2vec)))
        plt.plot(smooth(t2vec, smoothy), range(len(t2vec)), color='red', linewidth=3)
        plt.title('T2vec Plot')
        plt.xlabel('Value')
        plt.ylabel('Index')

        plt.tight_layout()
        plt.show()
    # Lines 9921-9924: Smooth vectors
    T1vecR = T1vec.copy()
    T2vecR = T2vec.copy()
    T1vec = smooth(T1vec, smoothy, edge_mirror=True)
    T2vec = smooth(T2vec, smoothy, edge_mirror=True)

    # Lines 9927-9941: Calculate parameters for first chunk
    if cur_junk == 1:
        # Calculate m1, m2, a, t0, maxis
        m1 = m2t(t1vec[0], a, t0, 0.5, SampInt, time=True)
        m2 = m2t(t2vec[0], a, t0, 0.5, SampInt, time=True)
        ai = (t2vec[0] * SampInt / 1e-10 - t1vec[0] * SampInt / 1e-10) / (np.sqrt(m2) - np.sqrt(m1))
        t0i = t1vec[0] * SampInt / 1e-10 - np.sqrt(m1) * ai
        
        if developing == 1:
            print([m1, m2, ai, t0i])

        m1 = m2t(t1vec[0], ai, t0i, 0.5, SampInt, time=True)
        m2 = m2t(t2vec[0], ai, t0i, 0.5, SampInt, time=True)
        a = (t2vec[0] * SampInt / 1e-10 - t1vec[0] * SampInt / 1e-10) / (np.sqrt(m2) - np.sqrt(m1))
        t0 = t1vec[0] * SampInt / 1e-10 - np.sqrt(m1) * ai
        
        if developing == 1:
            print([m1, m2, ai, t0i])
        
        maxis = m2t(np.arange(dims[0]), a, t0, 0.5, SampInt, time=True)


    # Lines 9950-9952: Set start and end indices
    ssttrrtt = 0 if cur_junk == 1 else smoothyhalf
    eenndd = len(T1vec) - 1 if cur_junk == n_junks else len(T1vec) - 2 - smoothyhalf

    # Lines 9955-9971: Main correction loop
    sumspecCORR = None
    avec = []
    t0vec = []
    for i in range(ssttrrtt, eenndd + 1):
        # Calculate ai, t0i, maxisi
        ai = (t2vec[i] * SampInt / 1e-10 - t1vec[i] * SampInt / 1e-10) / (np.sqrt(m2) - np.sqrt(m1))
        t0i = t1vec[i] * SampInt / 1e-10 - np.sqrt(m1) * ai
        avec.append(ai)
        t0vec.append(t0i)
        maxisi = m2t(np.arange(dims[0]), ai, t0i, 0.5, SampInt, time=True)

        # Interpolate and correct data
        sumspecCORRi = np.interp(maxis, maxisi, data[:, i])
        data[:, i] = sumspecCORRi

        # Accumulate sumspecCORR
        if i == ssttrrtt:
            sumspecCORR = sumspecCORRi
        else:
            sumspecCORR += sumspecCORRi

        if developing == 1:
            if i == smoothyhalf:
                plt.plot(sumspecCORRi[25000:30000])
            else:
                plt.plot(sumspecCORRi[25000:30000], color='red')
            print(np.ptp(sumspecCORRi))


    # Lines 9973-9976: Adjust data length
    if cur_junk == 1:
        length = len(T1vec) - smoothyhalf - 1
    elif cur_junk == n_junks:
        length = len(T1vec) - smoothyhalf
    else:
        length = len(T1vec) - smoothy
    ssttrrtt = 0 if cur_junk == 1 else 4

    data = data[:, ssttrrtt:ssttrrtt+length]

    # Lines 9979-9980: Create and return result structure
    s1 = {
        'sumspec': sumspecCORR,
        'data': data,
        'maxis': maxis,
        'avec': avec,
        't0vec': t0vec,
        'rest': rest,
        't1vec': T1vecR,
        't2vec': T2vecR,
        't1start': t1start,
        't2start': t2start,
        'm1': m1,
        'm2': m2
    }
    return s1

def masslib(extended=False, sulfate=False, cal=False):
    # Lines 10077-10084: Set up configuration based on input parameters
    ooo = 16 if extended else 5
    nnn = 2 if extended else 0
    Clll = 0
    sss = 1 if sulfate else 0
    if cal:
        SRI = 1
        ooo = 5
        nnn = 0
        Clll = 0
        sss = 0
    else:
        SRI = 0

    # Lines 10085-10103: Define atomic masses and constants
    C12, C13 = 12, 13.003355
    H1 = 1.007825
    O16, O17, O18 = 15.994915, 16.999131, 17.99916
    N14, N15 = 14.003074, 15.000108
    H_pos = 1.007276467
    e_neg = -0.000548533
    Cl35, Cl37 = 34.968852, 36.965903
    F = 18.998403
    S32O4 = 31.97207 + 4 * O16
    S34O4 = 33.967866 + 4 * O16
    Si28, Si29, Si30 = 27.976927, 28.976495, 29.97377
    I = 126.904473

    # Lines 10107-10314: Create library of protonated ion masses
    lib = []

    # Helper function to calculate mass and create entry
    def create_entry(n):
        mass = (C12*n[0] + C13*n[1] + H1*n[2] + O16*n[3] + O17*n[4] + O18*n[5] +
                N14*n[6] + N15*n[7] + Cl35*n[10] + Cl37*n[11] + F*n[12] +
                S32O4*n[13] + S34O4*n[14] + Si28*n[15] + Si29*n[16] + Si30*n[17] +
                I*n[18] + H_pos*n[8] + e_neg*n[9])
        return [mass] + n

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
    [1, 0, 0, 0, 0, 0, 0, 0, 0, 1, 3, 0, 0, 0, 0, 0, 0, 0, 0],  # CCl3+
    [1, 0, 0, 0, 0, 0, 0, 0, 0, 1, 2, 1, 0, 0, 0, 0, 0, 0, 0],  # CCl3+ isotope
    [0, 1, 0, 0, 0, 0, 0, 0, 0, 1, 3, 0, 0, 0, 0, 0, 0, 0, 0],  # CCl3+ isotope
    [0, 1, 0, 0, 0, 0, 0, 0, 0, 1, 2, 1, 0, 0, 0, 0, 0, 0, 0],  # CCl3+ isotope
    [1, 0, 0, 0, 0, 0, 0, 0, 0, 1, 1, 2, 0, 0, 0, 0, 0, 0, 0],  # CCl3+ isotope
    [6, 0, 5, 0, 0, 0, 0, 0, 1, 0, 0, 0, 0, 0, 0, 0, 0, 0, 1],  # Iodobenzene
    [5, 1, 5, 0, 0, 0, 0, 0, 1, 0, 0, 0, 0, 0, 0, 0, 0, 0, 1],  # Iodobenzene isotope
    [6, 0, 4, 0, 0, 0, 0, 0, 1, 0, 0, 0, 0, 0, 0, 0, 0, 0, 2],  # Di-Iodobenzene
    [5, 1, 4, 0, 0, 0, 0, 0, 1, 0, 0, 0, 0, 0, 0, 0, 0, 0, 2],  # Di-Iodobenzene isotope
    [6, 0, 4, 0, 0, 0, 0, 0, 1, 0, 0, 0, 0, 0, 0, 0, 0, 0, 1],  # Di-Iodobenzene fragment
    [5, 1, 4, 0, 0, 0, 0, 0, 1, 0, 0, 0, 0, 0, 0, 0, 0, 0, 1],  # Di-Iodobenzene fragment isotope
    [6, 0, 3, 0, 0, 0, 0, 0, 1, 0, 0, 0, 3, 0, 0, 0, 0, 0, 0],  # 3F-benzene
    [5, 1, 3, 0, 0, 0, 0, 0, 1, 0, 0, 0, 3, 0, 0, 0, 0, 0, 0],  # 3F-benzene isotope
    [6, 0, 3, 0, 0, 0, 0, 0, 1, 0, 3, 0, 0, 0, 0, 0, 0, 0, 0],  # 3Cl-benzene
    [6, 0, 3, 0, 0, 0, 0, 0, 1, 0, 2, 1, 0, 0, 0, 0, 0, 0, 0],  # 3Cl-benzene isotope
    [5, 1, 3, 0, 0, 0, 0, 0, 1, 0, 3, 0, 0, 0, 0, 0, 0, 0, 0],  # 3Cl-benzene isotope
    [5, 1, 3, 0, 0, 0, 0, 0, 1, 0, 2, 1, 0, 0, 0, 0, 0, 0, 0],  # 3Cl-benzene isotope
    [6, 0, 3, 0, 0, 0, 0, 0, 1, 0, 1, 2, 0, 0, 0, 0, 0, 0, 0],  # 3Cl-benzene isotope
]

    for ion in primary_ions:
        lib.append(create_entry(ion))

    # Lines 10186-10226: Add D and L ions
    for ss in range(4):
    # D ions
        lib.append(create_entry([6+ss*2, 0, 18+ss*6, 3+ss, 0, 0, 0, 0, 1, 0, 0, 0, 0, 0, 0, 3+ss, 0, 0, 0]))
        lib.append(create_entry([6+ss*2, 0, 18+ss*6, 3+ss, 0, 0, 0, 0, 1, 0, 0, 0, 0, 0, 0, 2+ss, 1, 0, 0]))
        lib.append(create_entry([6+ss*2, 0, 18+ss*6, 3+ss, 0, 0, 0, 0, 1, 0, 0, 0, 0, 0, 0, 2+ss, 0, 1, 0]))
        # D fragment -CH4
        lib.append(create_entry([5+ss*2, 0, 14+ss*6, 3+ss, 0, 0, 0, 0, 1, 0, 0, 0, 0, 0, 0, 3+ss, 0, 0, 0]))
        lib.append(create_entry([5+ss*2, 0, 14+ss*6, 3+ss, 0, 0, 0, 0, 1, 0, 0, 0, 0, 0, 0, 2+ss, 1, 0, 0]))
        lib.append(create_entry([5+ss*2, 0, 14+ss*6, 3+ss, 0, 0, 0, 0, 1, 0, 0, 0, 0, 0, 0, 2+ss, 0, 1, 0]))
        # L ions
        lib.append(create_entry([8+ss*2, 0, 24+ss*6, 2+ss, 0, 0, 0, 0, 1, 0, 0, 0, 0, 0, 0, 3+ss, 0, 0, 0]))
        lib.append(create_entry([8+ss*2, 0, 24+ss*6, 2+ss, 0, 0, 0, 0, 1, 0, 0, 0, 0, 0, 0, 2+ss, 1, 0, 0]))
        lib.append(create_entry([8+ss*2, 0, 24+ss*6, 2+ss, 0, 0, 0, 0, 1, 0, 0, 0, 0, 0, 0, 2+ss, 0, 1, 0]))
        # L fragment -CH4
        lib.append(create_entry([7+ss*2, 0, 20+ss*6, 2+ss, 0, 0, 0, 0, 1, 0, 0, 0, 0, 0, 0, 3+ss, 0, 0, 0]))
        lib.append(create_entry([7+ss*2, 0, 20+ss*6, 2+ss, 0, 0, 0, 0, 1, 0, 0, 0, 0, 0, 0, 2+ss, 1, 0, 0]))
        lib.append(create_entry([7+ss*2, 0, 20+ss*6, 2+ss, 0, 0, 0, 0, 1, 0, 0, 0, 0, 0, 0, 2+ss, 0, 1, 0]))




    # Lines 10232-10245: Add all other CHNOCl combinations
    for c in range(1, 41):
        for h in range(max(c-9, 0), c+2):
            for o in range(ooo+1):
                for nn in range(nnn+1):
                    for Cll in range(Clll+1):
                        for ss in range(sss+1):
                            if h > -0.5 and (nn >= 1 or h > 0.5):
                                n = [c, 0, 2*h+nn, o, 0, 0, nn, 0, 1, 0, Cll, 0, 0, ss, 0, 0, 0, 0, 0]
                                lib.append(create_entry(n))
                                if cal:  # charge transfer (SRI-mode)
                                    n[8], n[9] = 0, 1
                                    lib.append(create_entry(n))
                                if extended:  # 13C isotope
                                    n = [c-1, 1, 2*h+nn, o, 0, 0, nn, 0, 1, 0, Cll, 0, 0, ss, 0, 0, 0, 0, 0]
                                    lib.append(create_entry(n))

    # Line 10248: Sort the library by mass
    lib.sort(key=lambda x: x[0])

    # Line 10249: Print library size
    print(f"Library size: {len(lib)}")

    return lib

def masslib_old(extended=False, sulfate=False, cal=False):
    # Lines 10319-10329: Set configuration based on input parameters
    ooo = 16 if extended else 5
    nnn = 2 if extended else 0
    Clll = 0
    sss = 1 if sulfate else 0
    if cal:
        SRI = 1
        ooo = 5
        nnn = 0
        Clll = 0
        sss = 0
    else:
        SRI = 0

    # Lines 10330-10346: Define atomic masses and constants
    C12, C13 = 12, 13.003355
    H1 = 1.007825
    O16, O17, O18 = 15.994915, 16.999131, 17.99916
    N14, N15 = 14.003074, 15.000108
    H_pos = 1.007276467
    e_neg = -0.000548533
    Cl35, Cl37 = 34.968852, 36.965903
    F = 18.998403
    S32O4 = 31.97207 + 4 * O16
    S34O4 = 33.967866 + 4 * O16

    # Function to calculate mass from atom counts
    def calc_mass(n):
        return (C12*n[0] + C13*n[1] + H1*n[2] + O16*n[3] + O17*n[4] + O18*n[5] +
                N14*n[6] + N15*n[7] + Cl35*n[10] + Cl37*n[11] + F*n[12] +
                S32O4*n[13] + S34O4*n[14] + H_pos*n[8] + e_neg*n[9])

    # Lines 10350-10465: Define library of ions
    lib = []
    
    # Primary/inorganic ions
    ions = [
        [0, 0, 2, 0, 0, 1, 0, 0, 1, 0, 0, 0, 0, 0, 0],  # m21
        [0, 0, 4, 1, 0, 1, 0, 0, 1, 0, 0, 0, 0, 0, 0],  # m39
        [0, 0, 0, 2, 0, 0, 0, 0, 0, 1, 0, 0, 0, 0, 0],  # O2+
        [0, 0, 0, 1, 0, 1, 0, 0, 0, 1, 0, 0, 0, 0, 0],  # 16O18O+
        [0, 0, 0, 1, 0, 0, 1, 0, 0, 1, 0, 0, 0, 0, 0],  # NO+
        [0, 0, 0, 0, 0, 0, 2, 0, 1, 0, 0, 0, 0, 0, 0],  # N2H+
        [0, 0, 3, 0, 0, 0, 1, 0, 1, 0, 0, 0, 0, 0, 0],  # NH4+
        [0, 0, 6, 0, 0, 0, 2, 0, 1, 0, 0, 0, 0, 0, 0],  # NH3NH4+
        [0, 0, 1, 3, 0, 0, 1, 0, 1, 0, 0, 0, 0, 0, 0],  # HNO3H+
        [0, 0, 0, 2, 0, 0, 1, 0, 0, 1, 0, 0, 0, 0, 0],  # NO2+
        [0, 0, 0, 2, 0, 0, 0, 1, 0, 1, 0, 0, 0, 0, 0],  # 15NO2+
        [0, 0, 0, 1, 0, 1, 1, 0, 0, 1, 0, 0, 0, 0, 0],  # NO18O+
    ]
    
    for n in ions:
        lib.append([calc_mass(n), n])

    # CCl3+ and isotopes
    for n in [[1,0,0,0,0,0,0,0,0,1,3,0,0,0,0], [1,0,0,0,0,0,0,0,0,1,2,1,0,0,0],
              [0,1,0,0,0,0,0,0,0,1,3,0,0,0,0], [0,1,0,0,0,0,0,0,0,1,2,1,0,0,0],
              [1,0,0,0,0,0,0,0,0,1,1,2,0,0,0]]:
        lib.append([calc_mass(n), n])

    # C20H15ClO4 and isotopes
    for n in [[20,0,15,4,0,0,0,0,1,0,1,0,0,0,0], [19,1,15,4,0,0,0,0,1,0,1,0,0,0,0],
              [20,0,15,4,0,0,0,0,1,0,0,1,0,0,0], [19,1,15,4,0,0,0,0,1,0,0,1,0,0,0]]:
        lib.append([calc_mass(n), n])

    # C6H5ClO3 and isotopes
    for n in [[6,0,5,3,0,0,0,0,1,0,1,0,0,0,0], [5,1,5,3,0,0,0,0,1,0,1,0,0,0,0],
              [6,0,5,3,0,0,0,0,1,0,0,1,0,0,0], [5,1,5,3,0,0,0,0,1,0,0,1,0,0,0]]:
        lib.append([calc_mass(n), n])

    # Lines 10466-10472: Generate all other CHNOCl combinations
    for c in range(1, 41):
        for h in range(c-9, c+2):
            for o in range(ooo+1):
                for nn in range(nnn+1):
                    for Cll in range(Clll+1):
                        for ss in range(sss+1):
                            if h > -0.5 and (nn >= 1 or h > 0.5):
                                n = [c, 0, 2*h+nn, o, 0, 0, nn, 0, 1, 0, Cll, 0, 0, ss, 0]
                                lib.append([calc_mass(n), n])
                                if cal:  # charge transfer (SRI-mode)
                                    n[8], n[9] = 0, 1
                                    lib.append([calc_mass(n), n])
                                if extended:  # 13C isotope
                                    n = [c-1, 1, 2*h+nn, o, 0, 0, nn, 0, 1, 0, Cll, 0, 0, ss, 0]
                                    lib.append([calc_mass(n), n])

    # Sort the library by mass
    lib.sort(key=lambda x: x[0])

    print(len(lib))
    return lib

def resimodel(xmin, xmax, sumspec, baseline, resolution, resolution2nd, peakshape, a, t0, ex, sampint, instrument):
    # Line 10476: Smooth the sum spectrum
    smoo = SmoothSumSpec(sumspec, a, t0, ex, sampint, instrument)
    data = smoo.data
    
    # Lines 10478-10484: Prepare data and scales
    x = range(len(sumspec) - 1)
    xmin2 = xmin - (xmax - xmin) // 2
    xmax2 = xmax + (xmax - xmin) // 2
    baseline2 = baseline[xmin2:xmax2]
    data2 = data[xmin2:xmax2]
    mscale2 = m2t(x[xmin2:xmax2], a, t0, ex, sampint, time=True)

    # Lines 10486-10490: Calculate LOD5
    helft = (xmax2 - xmin2) // 2
    dataFIRSThalf = data2[:helft]
    dataSECONDhalf = data2[helft:]
    LOD5_1 = 3 * 2 * np.std(dataFIRSThalf[dataFIRSThalf < np.quantile(dataFIRSThalf, 0.25)])
    LOD5_2 = 3 * 2 * np.std(dataSECONDhalf[dataSECONDhalf < np.quantile(dataSECONDhalf, 0.25)])
    LOD5 = max(LOD5_1, LOD5_2)

    # Lines 10491-10493: Slice data and scales
    baseline = baseline[xmin:xmax]
    data = data[xmin:xmax]
    mscale = m2t(x[xmin:xmax], a, t0, ex, sampint, time=True)

    # Lines 10496-10584: Main processing loop
    length = len(data)
    modelCRUDE = np.zeros((3, 1))
    resi = data - baseline
    resi = data
    clean = 0

    for i in range(1, 13):
        # Peak detection and analysis
        PEAK = np.sum([resi, np.pad(resi[:-2], (0, 2)), np.pad(resi[1:-1], (1, 1)), 
                       np.pad(resi[:-3], (0, 3)), np.pad(resi[2:], (2, 0))], axis=0) / 5
        PEAK = np.argmax(PEAK)
        masss = mscale[PEAK]
        FWHM = masss / resolution
        SIGM = FWHM / (2 * np.sqrt(2 * np.log(2)))
        dm = max(0.002, sigm / 2)

        # Peakshape adjustment and interpolation
        sh1 = peakshape.copy()
        sh1[:, 0] = sh1[:, 0] * FWHM + masss
        sp1 = np.interp(mscale, sh1[:, 0], sh1[:, 1])

        # Calculate indices and factors
        indS = int((m2t(masss + sigm * 0.025, a, t0, ex, sampint, mass=True) - 
                    m2t(masss, a, t0, ex, sampint, mass=True)) / 2)
        indPts = 2 * indS + 1
        indd = range(PEAK - indS, PEAK + indS + 1)
        factor = np.sum(data[indd] - np.mean(baseline)) / np.sum(sp1[indd])
        factor = np.sum(resi[indd]) / np.sum(sp1[indd])

        # Model update and optimization
        modelCRUDE = np.column_stack((modelCRUDE, [masss, factor, resolution]))
        if i == 1:
            modelCRUDE = modelCRUDE[:, 1:]
            oldquality = 9999
        else:
            oldquality = S1.quality

        S1 = optimizer(data2, baseline2, mscale2, modelCRUDE, resolution, resolution2nd, peakshape)
        
        if i == 1:
            modelli = S1.model
        
        factors = S1.model
        masses = factors[0, :]
        masses.sort()
        mindiff = min(np.abs(np.diff(np.concatenate(([0], masses)))))
        factors = factors[1, :]

        if min(factors) < LOD5 or i == 12:
            clean = 1
        else:
            modelli = S1.model

        # Clean up and finalize model
        if clean == 1:
            hlppp = np.argsort(modelli[0, :])
            ff = len(hlppp)
            modelli2 = modelli[:, hlppp]  # m/z sorted
            modelli = modelli2[:, 0:1]
            
            if ff > 1:
                for qq in range(1, ff):
                    LastIndex = np.argmax(modelli[0, :])
                    if modelli2[0, qq] - np.max(modelli[0, :]) < dm:
                        modelli[1, LastIndex] += modelli2[1, qq]
                    elif modelli2[0, qq] - np.max(modelli[0, :]) > 2 * dm:
                        modelli = np.column_stack((modelli, modelli2[:, qq:qq+1]))
                    elif max(modelli[1, LastIndex], modelli2[1, qq]) / min(modelli[1, LastIndex], modelli2[1, qq]) < 5:
                        modelli = np.column_stack((modelli, modelli2[:, qq:qq+1]))
                    else:
                        if modelli[1, LastIndex] < modelli2[1, qq]:
                            modelli[0, LastIndex] = modelli2[0, qq]
                        modelli[1, LastIndex] += modelli2[1, qq]

            S1 = optimizer(data2, baseline2, mscale2, modelli, resolution, resolution2nd, peakshape)
            modelli = S1.model
            i = 100

        # Print results and update model
        print(f"{LOD5:7.1f}{1000*dm:7.1f}    residue: {S1.quality:7.2f} permil,    improvement: {100-100*S1.quality/oldquality:7.2f}%")
        hlppp = np.argsort(modelli[0, :])
        print(modelli[:, hlppp])
        print("  END   ")
        reconstruction = reconstruct(modelli, baseline, mscale, resolution, peakshape)
        resi = data - reconstruction[:, 0]
        modelCRUDE = modelli.copy()

    # Final reconstruction and model creation
    priori = reconstruct(modelli, baseline, mscale, resolution, peakshape)
    posteriori = reconstruct(modelli, baseline, mscale, resolution, peakshape)
    modelpeaks = {
        'mass': masss,
        'bin': m2t(masss, a, t0, ex, sampint, mass=True),
        'binscale': x[xmin:xmax],
        'mscale': mscale,
        'posteriori': posteriori,
        'residual': data - posteriori,
        'priori': priori,
        'model': modelli,
        'LOD5': LOD5,
        'dm': dm
    }

    return modelpeaks

def reconstruct(model, baseline, mscale, resolution, peakshape):
    # Line 10586-10588: Initialize reconstruction arrays
    reconstruction = baseline.copy()
    reconstruction2 = baseline.copy()
    reconstruction3 = baseline.copy()

    # Line 10589-10590: Determine the number of peaks to process
    n = model.shape[1] if len(model.shape) > 1 else 1

    # Line 10591-10605: Loop through each peak and add to reconstructions
    for i in range(n):
        # Extract peak parameters
        masss = model[0, i]
        factor = model[1, i]
        resolution3 = model[2, i]

        # Adjust peakshape for current mass and resolution
        sh1 = peakshape.copy()
        sh1[:, 0] = (sh1[:, 0] * masss / resolution3 + masss)

        # Interpolate peak shape components to match mscale
        sp1 = interp1d(sh1[:, 0], sh1[:, 1], kind='linear', fill_value='extrapolate')(mscale)
        sp2 = interp1d(sh1[:, 0], sh1[:, 2], kind='linear', fill_value='extrapolate')(mscale)
        sp3 = interp1d(sh1[:, 0], sh1[:, 3], kind='linear', fill_value='extrapolate')(mscale)

        # Add scaled peak shapes to reconstructions
        reconstruction += sp1 * factor
        reconstruction2 += sp2 * factor
        reconstruction3 += sp3 * factor

    # Line 10606: Return the three reconstructions as a 2D array
    return np.column_stack((reconstruction, reconstruction2, reconstruction3))

def optimizer(data2, baseline2, mscale2, model, resolution, resolution2nd, peakshape):
    # Line 10608-10610: Initialize quality variables
    qualityold = 0
    qualityold2 = -1
    resmult = 0.99

    # Lines 10611-10651: Main optimization loop
    for iteration in range(501):
        reconstruction = reconstruct(model, baseline2, mscale2, resolution, peakshape)
        reconstruction = reconstruction[:, 0]
        quality = 1000.0 * sum(abs(data2 - reconstruction)) / sum(abs(data2))
        qualityold = [quality] + qualityold[:9]

        if iteration > 9:
            modelQ = ResidueQuality(data2, baseline2, reconstruction, model, mscale2, resolution)
            offset = modelQ[5, :]
            asym = modelQ[4, :] - modelQ[6, :]
            ioff = np.argmax(abs(offset))
            iasy = np.argmax(abs(asym))
            model[1, ioff] += offset[ioff] / 2
            model[1, ioff] = max(model[1, ioff], 0)
            model[0, iasy] -= model[0, iasy] * asym[iasy] / (resolution * max(8000, 8000 * abs(asym[iasy]) / 50))

        # Lines 10652-10688: Check for convergence and adjust resolution
        if iteration > 9 and (abs(quality - np.mean(qualityold)) < 0.005 or abs(np.mean(qualityold2) - np.mean(qualityold)) < 0.005):
            dims = modelQ.shape
            nPEAKS = 1 if len(dims) == 1 else dims[1]
            for sss in range(nPEAKS):
                if modelQ[4, sss] < modelQ[5, sss] and modelQ[6, sss] < modelQ[5, sss]:
                    resolution2 = resolution
                    for iteration2 in range(300):
                        resolution2 /= resmult
                        model[2, sss] = resolution2
                        reconstruction = reconstruct(model, baseline2, mscale2, resolution2, peakshape)
                        reconstruction = reconstruction[:, 0]
                        quality = 1000.0 * sum(abs(data2 - reconstruction)) / sum(abs(data2))
                        qualityold = [quality] + qualityold[:9]
                        modelQ = ResidueQuality(data2, baseline2, reconstruction, model, mscale2, resolution2)
                        print([resolution2] + list(modelQ[:, sss]))
                        if (modelQ[4, sss] > modelQ[5, sss] and modelQ[6, sss] > modelQ[5, sss]) or resolution2 > resolution * 1.15:
                            model[2, sss] = resolution if resolution2 > resolution * 1.15 else resolution2
                            break
                elif modelQ[4, sss] > modelQ[5, sss] and modelQ[6, sss] > modelQ[5, sss]:
                    resolution2 = resolution
                    for iteration2 in range(300):
                        resolution2 *= resmult
                        model[2, sss] = resolution2
                        reconstruction = reconstruct(model, baseline2, mscale2, resolution2, peakshape)
                        reconstruction = reconstruction[:, 0]
                        quality = 1000.0 * sum(abs(data2 - reconstruction)) / sum(abs(data2))
                        qualityold = [quality] + qualityold[:9]
                        modelQ = ResidueQuality(data2, baseline2, reconstruction, model, mscale2, resolution2)
                        print([resolution2] + list(modelQ[:, sss]))
                        if modelQ[4, sss] < modelQ[5, sss] or modelQ[6, sss] < modelQ[5, sss]:
                            model[2, sss] = resolution2 if resolution2 < resolution2nd + (resolution - resolution2nd) * 0.3 else resolution
                            break
            break

        # Lines 10689-10690: Print iteration results
        print([iteration, quality, np.mean(qualityold2), np.mean(qualityold)])
        qualityold2 = qualityold

    # Lines 10698-10700: Return optimized model and quality
    return {'model': model, 'quality': quality}


def ResidueQuality(data, baseline, reconstruction, model, mscale, resolution):
    # Line 10699: Determine the number of model components
    n = len(model)
    if isinstance(n, (list, tuple)):
        n = max(n[1])
    else:
        n = 1

    # Line 10701: Initialize modelQ array
    modelQ = [[0.0] * 7 for _ in range(n)]  # mass, factor, Max(data_center), Mean(data_center), Qleft, Qcenter, Qright

    # Line 10702: Copy model data to modelQ
    for i in range(n):
        modelQ[i][0:3] = model[i][0:3]

    # Lines 10704-10726: Main loop for quality calculations
    for i in range(n):
        masss = modelQ[i][0]
        factor = modelQ[i][1]
        resolution3 = modelQ[i][2]
        FWHM = masss / resolution3

        # Find indices for integration regions
        i1 = max([j for j, m in enumerate(mscale) if m < masss - FWHM/2])
        i2 = max([j for j, m in enumerate(mscale) if m < masss - FWHM/6])
        i3 = max([j for j, m in enumerate(mscale) if m < masss + FWHM/6])
        i4 = max([j for j, m in enumerate(mscale) if m < masss + FWHM/2])

        # Calculate quality metrics
        modelQ[i][2] = max(data[i2:i3] - baseline[i2:i3])
        modelQ[i][3] = sum(data[i2:i3] - baseline[i2:i3]) / (i3 - i2)
        modelQ[i][4] = (sum(data[i1:i2]) - sum(reconstruction[i1:i2])) / (i2 - i1)
        modelQ[i][5] = (sum(data[i2:i3]) - sum(reconstruction[i2:i3])) / (i3 - i2)
        modelQ[i][6] = (sum(data[i3:i4]) - sum(reconstruction[i3:i4])) / (i4 - i3)

    # Line 10731: Return the calculated quality metrics
    return modelQ

def model_rest(xmin, xmax, sumspec, baseline, resolution, peakshape, a, t0, ex, sampint, instrument):
    # Line 10735: Initialize factor2
    factor2 = factor.copy()

    # Lines 10736-10816: Main optimization loop
    for l in range(51):
        sensi = np.zeros((length2, 3))
        step = 0.00001
        modpeaks = np.zeros((len(mscale), length2, 3))
        Pcenter = mscale[PEAKS] + offset
        FWHM = Pcenter / (resolution * 1.0)
        SIGM = FWHM / (2 * np.sqrt(2 * np.log(2)))

        # Line 10743-10745: Update factor2
        if np.max(factor2) == 0:
            factor2 = factor.copy()
        jot = np.argmax(factor2)
        factor2[jot] = 0

        # Lines 10758-10784: Optimize offset
        fertig = 0
        counter = 0
        while fertig == 0:
            Pcenter = mscale[PEAKS] + offset
            FWHM = Pcenter / (resolution * 1.0)
            for i in range(length2):
                sh1 = peakshape.copy()
                sh1[:, 0] = (sh1[:, 0] * FWHM[i] + Pcenter[i])
                modpeaks[:, i, 1] = np.interp(mscale, sh1[:, 0], sh1[:, 1], left=np.nan, right=np.nan) * factor[i]

            model = np.sum(modpeaks, axis=1) + base if length2 > 1.5 else modpeaks.reshape(-1) + base
            residual = data - model[:, 1]

            if l == 0 and counter == 0:
                priori = residual.copy()

            if counter > 0.5:
                qualit_prev = qualit
            else:
                offset[jot] += step

            qualit = np.mean(np.abs(residual)) / np.median(np.abs(residual))

            if counter > 0.5:
                if qualit < qualit_prev:
                    offset[jot] += step
                elif counter < 1.5:
                    step = -step
                    offset[jot] += step
                else:
                    fertig = 1

            counter += 1

        print([Pcenter[jot], offset * 1000])

        # Lines 10787-10814: Optimize factor
        fertig = 0
        counter = 0
        fstep = 1.0001
        while fertig == 0:
            Pcenter = mscale[PEAKS] + offset
            FWHM = Pcenter / (resolution * 1.0)
            for i in range(length2):
                sh1 = peakshape.copy()
                sh1[:, 0] = (sh1[:, 0] * FWHM[i] + Pcenter[i])
                modpeaks[:, i, 1] = np.interp(mscale, sh1[:, 0], sh1[:, 1], left=np.nan, right=np.nan) * factor[i]

            model = np.sum(modpeaks, axis=1) if length2 > 1.5 else modpeaks.reshape(-1)
            residual = data - model[:, 1]

            if counter > 0.5:
                qualit_prev = qualit
            else:
                factor[jot] *= fstep

            qualit = np.mean(np.abs(residual)) / np.median(np.abs(residual))

            if counter > 0.5:
                if qualit < qualit_prev:
                    factor[jot] *= fstep
                elif counter < 1.5:
                    fstep = 1.0 / fstep
                    factor[jot] *= fstep
                else:
                    fertig = 1

            counter += 1

        print([Pcenter[jot], factor])

    # Lines 10817-10819: Create and return modelpeaks structure
    modelpeaks = {
        'mass': Pcenter,
        'bin': m2t(Pcenter, a, t0, ex, SampInt, mass=True),
        'binscale': binscale,
        'mscale': mscale,
        'model': model[:, 1],
        'residual': residual,
        'priori': priori
    }
    return modelpeaks

def overlap(masslist, data, name, path):
    # Line 10819: Start timing
    mistt = time.time()

    # Lines 10820-10827: Initialize variables and reshape data if necessary
    lenPeaks = len(masslist)
    dims = data.shape
    if dims[0] != lenPeaks:
        data = data.T
        lenCyc = dims[0]
    else:
        lenCyc = dims[1]
    
    if not isinstance(name, list):
        name = [name] * lenCyc
    
    nam = 'xxx'
    DataCorr = np.zeros((lenPeaks, lenCyc))

    # Lines 10829-10879: Main loop over cycles
    for j in range(lenCyc):
        if name[j] != nam:
            nam = name[j]
            # Read resolution and peak shape data
            reso = readfloat(f"{path}FileInfo/Time2MassPAR/{nam}Par2.fdt")[8]
            PeakData = peaktable(masslist, reso)
            FWHM = PeakData[:, 3] * 2 * np.sqrt(2 * np.log(2))
            Peakshape = readfloat(f"{path}FileInfo/PeakShape/{nam}PeakShape.fdt")
            x7 = Peakshape[:, 0]
            dataX3 = Peakshape[:, 1]

            # Calculate peak boundaries
            bb1 = (PeakData[:, 1] - PeakData[:, 0]) / FWHM
            bb2 = (PeakData[:, 2] - PeakData[:, 0]) / FWHM
            bb1_ = np.concatenate(([PeakData[:, 1], 0] - [0, PeakData[:, 0]]) / [FWHM, 1])
            bb2_ = (PeakData[:, 2] - np.concatenate([PeakData[1:, 0], [10000]])) / FWHM
            b1, b2, b1_, b2_ = bb1.copy(), bb2.copy(), bb1_.copy(), bb2_.copy()

            # Calculate peak integrals
            for i in range(len(b1)):
                b1[i] = np.sum(dataX3[x7 < bb1[i]]) / np.sum(dataX3)
                b2[i] = np.sum(dataX3[x7 > bb2[i]]) / np.sum(dataX3)
                filter_ = x7 > bb1_[i]
                b1_[i] = np.sum(dataX3[filter_]) / np.sum(dataX3) if filter_.any() else 0
                filter_ = x7 < bb2_[i]
                b2_[i] = np.sum(dataX3[filter_]) / np.sum(dataX3) if filter_.any() else 0

            # Build matrix A
            dim = len(b1)
            A = np.zeros((dim, dim))
            for i in range(dim):
                if i > 0:
                    A[i, i-1] = b2_[i-1]
                A[i, i] = 1 - b1[i] - b2[i]
                if i < dim - 1:
                    A[i, i+1] = b1_[i+1]

            print(f"Matrix A built, dims: {A.shape}, elapsed time: {time.time() - mistt}")

            # Prepare for LU decomposition
            C = A.copy()
            lu, piv = scipy.linalg.lu_factor(C)

        # Solve the system for each cycle
        B = data[:, j]
        B2 = np.zeros(dim)
        B2[B > -9999] = B[B > -9999]
        S = scipy.linalg.lu_solve((lu, piv), B2)
        DataCorr[:, j] = S

        if j % 200 == 0:
            print(f'Cycle: {j}, time = {time.time() - mistt}')

    # Lines 10881-10882: Transpose DataCorr if necessary
    if dims[0] != lenPeaks:
        DataCorr = DataCorr.T

    return DataCorr

def PS1(event, SumSpectrum, SampInt, duration, cycles, extractions, destfold, name, instrument):
    # Lines 10889-10891: Update widget info
    WIDGET_CONTROL(WIDGET_INFO(event.TOP, FIND_BY_UNAME='Text_info'), SET_VALUE='Peaks & time2mass:')
    info = WIDGET_CONTROL(WIDGET_INFO(event.TOP, FIND_BY_UNAME='Text_info'), GET_VALUE=info)

    # Lines 10895-10896: Initialize variables
    x = range(len(SumSpectrum) - 1)
    name = name.reshape(-1)
    lib = masslib(cal=True)[0, :]

    # Lines 10899-10900: Get ion parameters
    ions = [getpar('ION1'), getpar('ION2'), getpar('ION3')]
    devs = [getpar('DEV1'), getpar('DEV2'), getpar('DEV3')]

    # Lines 10903-10913: Set instrument-specific parameters
    if instrument.startswith('TOF1000'):
        a = getpar('default_a_1000')
        t0 = getpar('default_t0_1000')
    elif instrument.startswith('TOF8000'):
        a = getpar('default_a_8000')
        t0 = getpar('default_t0_8000')
    elif instrument.startswith('VOCUS'):
        a = getpar('default_a_8000')
        t0 = getpar('default_t0_8000')
    print([0, 0, 0, a, t0])

    # Line 10914: Detect peaks
    peaklist = DetectPeaks(SumSpectrum, a, t0, 0.5, SampInt, instrument)

    # Line 10917: Calibrate crude
    par = CalCrude(peaklist, SampInt, lib, instrument)

    # Lines 10920-10926: Use crude or default calibration
    if getpar('useCRUDEdefault') == 0:
        a = par.a
        t0 = par.t0
    else:
        a = getpar('CRUDEa')
        t0 = getpar('CRUDEt0')

    # Lines 10929-10933: Set exponent and update widget info
    ex = 0.5
    WIDGET_CONTROL(WIDGET_INFO(event.TOP, FIND_BY_UNAME='Text_info'), 
                   SET_VALUE=[info, 'CalCrude:', f'a={a}', f't0={t0}', 'ex=0.5'])

    # Line 10935: Get top 16 peaks
    Top16 = par.Top16

    # Lines 10936-11530: Main processing block
    if a != -9999:  # if calCrude didn't fail
        # Detect peaks with clean option
        peaklist = DetectPeaks(SumSpectrum, a, t0, ex, SampInt, instrument, clean=True)
        ex = 0.5
        a_raw, t0_raw, ex_raw = a, t0, 0.5

        # Calculate resolution and masses
        reso = np.mean(peaklist[:, 3])
        masses = m2t(peaklist[:, 7], a_raw, t0_raw, 0.5, SampInt, time=True)

        # Define mass ranges and calculate signals
        mm_1 = [21.022, 39.0327, 30.994, 33.9935, 47.9966]
        mm_2 = [487, 242, 269, 242, 242]
        Psig = [0]
        for r in range(5):
            # Calculate signal and background for each mass range
            rr2 = np.floor(m2t([mm_1[r] - max(0.005, mm_1[r]/reso), mm_1[r] + max(0.005, mm_1[r]/reso)], a, t0, ex, SampInt))
            sig2 = (corrtr(79)/corrtr(mm_1[r])) * mm_2[r] * np.sum(SumSpectrum[int(rr2[0]):int(rr2[1])]) / duration if min(rr2) > 0 else 0
            pp = int(np.floor((rr2[1] - rr2[0]) / 2))
            qq = pp - 1 if int(np.floor(1 + (rr2[1] - rr2[0]) / 2)) == pp else pp
            bk2 = (corrtr(79)/corrtr(mm_1[r])) * mm_2[r] * np.sum([SumSpectrum[int(-1-qq+rr2[0]):int(rr2[0]-1)], SumSpectrum[int(1+rr2[1]):int(rr2[1]+1+pp)]]) / duration if min(rr2) > pp and sig2 > 0 else 0
            Psig.append(sig2 - bk2)
        Psig = Psig[1:]

        # Determine ionization mode
        mode = 0 if np.argmax(Psig) < 1.5 else (9 if np.argmax(Psig) == 3 else 5)

        # Filter library based on mode
        Flib = -1
        if mode == 0:
            Flib = np.where((2 * np.floor(np.floor(lib + 0.4) / 2) != np.floor(lib + 0.4)) | (lib < 50))[0]
            lib2 = lib[Flib] if max(Flib) > 0 else lib
        else:
            lib2 = lib

        # Perform fine calibration
        fineCal = CalFine(peaklist[:, 7], lib2, a, t0, SampInt, mode, instrument)
        a_fine, t0_fine, ex_fine = max(fineCal.a), max(fineCal.t0), max(fineCal.ex)
        rat120 = fineCal.rat120

        # Update widget info
        WIDGET_CONTROL(WIDGET_INFO(event.TOP, FIND_BY_UNAME='Text_info'), GET_VALUE=info)
        WIDGET_CONTROL(WIDGET_INFO(event.TOP, FIND_BY_UNAME='Text_info'), 
                       SET_VALUE=[info, '', 'CalFine:', f'a={a_fine}', f't0={t0_fine}', f'ex={ex_fine}'])

        # Perform 3-point calibration
        cal3 = Cal3pt(peaklist, a_fine, t0_fine, ex_fine, SampInt, lib2, mode, instrument)
        a3, t03, ex3 = max(cal3.a), max(cal3.t0), max(cal3.ex)

        # Update widget info
        WIDGET_CONTROL(WIDGET_INFO(event.TOP, FIND_BY_UNAME='Text_info'), GET_VALUE=info)
        WIDGET_CONTROL(WIDGET_INFO(event.TOP, FIND_BY_UNAME='Text_info'), 
                       SET_VALUE=[info, '', 'Cal3pt:', f'a={a3}', f't0={t03}', f'ex={ex3}'])

        # Test scale for different calibrations
        tst = testSCALE(m2t(peaklist[:, 7], a_raw, t0_raw, 0.5, SampInt, time=True), lib2, instrument, ions, devs)
        scorRaw = tst.scorppm
        tst = testSCALE(m2t(peaklist[:, 7], a_fine, t0_fine, ex_fine, SampInt, time=True), lib2, instrument, ions, devs)
        scorParVar = tst.scorppm
        tst = testSCALE(m2t(peaklist[:, 7], a3, t03, ex3, SampInt, time=True), lib2, instrument, ions, devs)
        scor3 = tst.scorppm

        # Select best calibration
        if scorRaw > scorParVar and scorRaw > scor3:
            a, t0, ex = a_raw, t0_raw, ex_raw
        elif scor3 > scorParVar and scor3 > scorRaw:
            a, t0, ex = a3, t03, ex3
        elif scorParVar >= scor3 and scorParVar >= scorRaw:
            a, t0, ex = a_fine, t0_fine, ex_fine

        if getpar('ForceThroughMaSet') == 1:
            a, t0, ex = a3, t03, ex3
        # Baseline signal calculation
        if getpar('PeakAnalysis') == 1:
            times = peaklist[:, 7]
            masses = m2t(times, a, t0, ex, SampInt, time=True)
            extr = float(cycles) * float(extractions)
            SumSpectrumSM = smooth(SumSpectrum, 5 * floor(0.55 + 1e-10 / SampInt))
            endmass = m2t(len(SumSpectrum), a, t0, ex, SampInt)

            # Calculate baseline signal
            pts = len(SumSpectrum)
            baseline = np.zeros(pts)
            oldt = 0
            n_SEG = floor(600 * 1e-10 / SampInt)
            
            for iii in range(floor(pts / n_SEG)):
                newt = (iii + 1) * n_SEG
                SEG = SumSpectrum[oldt:newt]
                for iv in range(7):
                    nn = len(SEG)
                    i_max = np.argmax(SEG)
                    if i_max <= n_SEG // 20:
                        SEG = SEG[n_SEG // 10:]
                    elif i_max >= nn - n_SEG // 20:
                        SEG = SEG[:-(n_SEG // 10)]
                    else:
                        SEG = np.concatenate((SEG[:i_max - n_SEG // 20], SEG[i_max + n_SEG // 20:]))
                baseline[oldt:newt] = np.mean(SEG)
                oldt = newt
            
            baseline[baseline == 0] = np.min(baseline[baseline > 0])
            baselineSM = np.convolve(baseline, savgol_filter(floor(1.5 * n_SEG), floor(1.5 * n_SEG), 0, 1), mode='same')

            # Plotting
            plt.figure(figsize=(11.5, 7))
            plt.subplot(411)
            plot_baseline(SumSpectrum, baselineSM, extr, 2e4, 14e4)
            plt.subplot(412)
            plot_baseline(SumSpectrum, baselineSM, extr, 14e4, 22e4)
            plt.subplot(413)
            plot_baseline(SumSpectrum, baselineSM, extr, 22e4, 28e4)
            plt.subplot(414)
            plot_baseline(SumSpectrum, baselineSM, extr, 28e4, pts - 1)
            plt.tight_layout()
            plt.savefig(f'{destfold}FileInfo/Baseline/Baseline_{name[0]}.jpg')
            plt.close()

        def plot_baseline(SumSpectrum, baselineSM, extr, start, end):
            x = np.arange(start, end)
            y = SumSpectrum[start:end] / extr
            plt.semilogy(x, y, 'k-')
            plt.semilogy(x, baselineSM[start:end] / extr, 'r-', linewidth=2)
            plt.ylim(np.min(baselineSM / extr), 5 * np.max(baselineSM[start:end] / extr))

        
        # Output preparation
        maxMass = m2t(max(x), a, t0, ex, SampInt)
        list = [18.0338, 21.0221, 29.0134, 29.9974, 31.9893, 38.0326, 59.0491, 79.0542, 93.0699, 137.1325]
        masses = m2t(peaklist[:, 7], a, t0, ex, SampInt, time=True)
        res = [0]
        counts = [0]
        for k in range(10):
            cnts = 0
            rs = 0
            if min(abs(masses - list[k])) < 0.003:
                i1 = np.argmin(abs(masses - list[k]))
                tL = m2t(m2t(peaklist[i1, 7], a, t0, ex, SampInt, time=True) * (1 - 1/(reso*2*np.sqrt(2*np.log(2)))), a, t0, ex, SampInt)
                tH = m2t(m2t(peaklist[i1, 7], a, t0, ex, SampInt, time=True) * (1 + 1/(reso*2*np.sqrt(2*np.log(2)))), a, t0, ex, SampInt)
                subdata = SumSpectrum[int(tL)-1:int(tH)+2]
                subx = x[int(tL)-1:int(tH)+2]
                coeff = gaussfit(subx, subdata, nterms=3)
                cnts = coeff[0] * coeff[2] * np.sqrt(2*np.pi)
                ms = coeff[1]
                msL = coeff[1] - np.sqrt(2*np.log(2)) * coeff[2]
                msH = coeff[1] + np.sqrt(2*np.log(2)) * coeff[2]
                rs = ((ms-t0)/a)**2 / (((msH-t0)/a)**2 - ((msL-t0)/a)**2)
            else:
                cnts = -9999
                rs = -9999
            counts.append(cnts)
            res.append(rs)

        output = [maxMass, a, t0] + counts[1:11] + res[1:11]
        massnames = ['m18.0338', 'm21.0221', 'm29.0134', 'm29.9974', 'm31.9893', 'm38.0326', 'm59.0491', 'm79.0542', 'm93.0699', 'm137.1325']
        peaklist[:, 0] = masses
        peaklist[:, 1] = m2t(peaklist[:, 1], a, t0, ex, SampInt, time=True)
        peaklist[:, 2] = m2t(peaklist[:, 2], a, t0, ex, SampInt, time=True)
        peaklist[:, 3] = masses / (peaklist[:, 2] - peaklist[:, 1])
        peaklistnames = ['PeakMax[Da]', 'StartPeak[Da]', 'EndPeak[Da]', 'PeakBroadness', 'SlopeMax', 'SlopeMin', 'PeakHight[counts]', 'PeakMax[time]', 'BaseHight[counts]']

        # Create output structure
        s1 = {
            'MaxMass': maxMass, 'a': a, 't0': t0, 'ex': ex, 'scorParVar': scorParVar,
            'a3': a3, 't03': t03, 'ex3': ex3, 'scor3': scor3, 'massnames': massnames,
            'counts': counts[1:11], 'res': res[1:11], 'masslistnames': peaklistnames,
            'masslist': peaklist, 'baseline': baselineSM, 
            'PeakShape': np.column_stack((xPeakshape, Peakshape, Peakshape2, Peakshape3)),
            'resolution': reso, 'mode': mode
        }
    else:
        # If calCrude failed, create a structure with default values
        s1 = {
            'MaxMass': -9999, 'a': -9999, 't0': -9999, 'ex': -9999, 'scorParVar': -9999,
            'a3': -9999, 't03': -9999, 'ex3': -9999, 'scor3': -9999, 'massnames': -9999,
            'counts': -9999, 'res': -9999, 'masslistnames': -9999, 'masslist': -9999,
            'baseline': -9999, 'PeakShape': -9999, 'resolution': -9999, 'mode': -9999
        }

    return s1

def PS2(x, y):
    # Line 11535: Set separation limit
    sepp = 8  # limit to separate independent peaks, unit: bins of bin mass scale, i.e. 64 ppm above 125 Da
    
    # Lines 11536-11539: Smooth data and calculate derivatives
    ysm = smooth(y, 5)
    dysm = deriv(x, ysm)
    length = len(y)
    index = range(length)
    
    # Lines 11540-11545: Find peaks
    PEAKS = np.maximum.reduce([ysm, np.pad(ysm[:-1], (1,0)), np.pad(ysm[1:], (0,1)),
                               np.pad(ysm[:-2], (2,0)), np.pad(ysm[2:], (0,2))])
    PEAKS[dysm < 0] = 0
    PEAKS[np.pad(dysm[1:], (0,1)) > 0] = 0
    PEAKS[PEAKS < max(ysm) * getpar('LowThresh')] = 0
    PEAKS[PEAKS < 0.55] = 0
    PEAKS = [i for i in index if PEAKS[i] > 0]
    length2 = len(PEAKS)  # number of peaks
    
    # Lines 11547-11556: Initialize arrays
    print('HHGHGH', length2)
    xpeaks, sgma, ymax, ymax2, aerM1s, aerM2s, aerP1s, aerP2s = [0], [0], [0], [0], [0], [0], [0], [0]
    notres = [1]
    
    # Lines 11557-11610: Main loop for peak analysis
    for i in range(length2):
        indis = range(PEAKS[i] - 5, PEAKS[i] + 6)
        xi = [x[j] for j in indis]
        
        # Gaussian fitting and parameter calculation
        if max([y[j] for j in indis]) <= 10:
            coeff = gaussfit(xi, [y[j] for j in indis], nterms=3)
        else:
            coeff = gaussfit(xi, [y[j] + 0.0001 for j in indis])
        
        # Additional fitting if necessary
        if coeff[0] < 0 or coeff[0] > max([y[j] for j in indis]) * 3 or 1e6 * abs(coeff[2]) / mean(xi) < 5 or coeff[2] <= 0:
            coeff = gaussfit(xi, [ysm[j] for j in indis], nterms=3)
        
        # Calculate separation limit and update arrays
        seplim = max(float(sepp) * 1e6 * (x[next(j for j in range(len(x)) if abs(x[j] - coeff[1]) == min(abs(x[k] - coeff[1]) for k in range(len(x))))] - 
                                          x[next(j for j in range(len(x)) if abs(x[j] - coeff[1]) == min(abs(x[k] - coeff[1]) for k in range(len(x))))]) / coeff[1])
        xpeaks.append(coeff[1])
        ymax.append(coeff[0])
        flt = next(j for j in range(len(x)) if abs(x[j] - coeff[1]) == min(abs(x[k] - coeff[1]) for k in range(len(x))))
        ymax2.append(ysm[flt])
        sgma.append(coeff[2])
        
        # Calculate areas
        flt2 = [j for j in range(len(x)) if coeff[1] - 2 * coeff[2] < x[j] <= coeff[1] - coeff[2]]
        aerM2s.append(sum([y[j] for j in flt2]))
        flt2 = [j for j in range(len(x)) if coeff[1] - coeff[2] < x[j] <= coeff[1]]
        aerM1s.append(sum([y[j] for j in flt2]))
        flt2 = [j for j in range(len(x)) if coeff[1] < x[j] <= coeff[1] + coeff[2]]
        aerP1s.append(sum([y[j] for j in flt2]))
        flt2 = [j for j in range(len(x)) if coeff[1] + coeff[2] < x[j] <= coeff[1] + 2 * coeff[2]]
        aerP2s.append(sum([y[j] for j in flt2]))
        
        # Determine if peak should be kept
        if coeff[0] < 0 or coeff[0] > 3.5 * max([y[j] for j in indis]):
            notres.append(1)
        elif coeff[1] < min(xi) or coeff[1] > max(xi):
            notres.append(1)
        elif coeff[2] < 0:
            notres.append(1)
        elif min(abs(xp - coeff[1]) / coeff[1] for xp in xpeaks[:i+1]) < seplim * 1e-6:
            check = 1
        else:
            notres.append(0)
        
        # Handle overlapping peaks
        if check == 1:
            iii = [j for j in range(i+1) if 1e6 * abs((xpeaks[j] - coeff[1]) / coeff[1]) < seplim]
            if min(notres[j] for j in iii) == 1:
                notres.append(0)
            else:
                iv = [j for j in iii if notres[j] == 0]
                if max(ymax2[j] for j in iv) > ysm[flt]:
                    notres.append(1)
                else:
                    for j in iv:
                        notres[j] = 1
                    notres.append(0)
    
    # Lines 11612-11622: Post-processing
    xpeaks = [xpeaks[i] for i in range(len(notres)) if notres[i] == 0]
    sgma = [sgma[i] for i in range(len(notres)) if notres[i] == 0]
    ymax = [ymax[i] for i in range(len(notres)) if notres[i] == 0]
    ymax2 = [ymax2[i] for i in range(len(notres)) if notres[i] == 0]
    aerM1s = [aerM1s[i] for i in range(len(notres)) if notres[i] == 0]
    aerM2s = [aerM2s[i] for i in range(len(notres)) if notres[i] == 0]
    aerP1s = [aerM1s[i] for i in range(len(notres)) if notres[i] == 0]
    aerP2s = [aerM2s[i] for i in range(len(notres)) if notres[i] == 0]
    PEAKS = [0] + [PEAKS[i] for i in range(len(notres)) if notres[i] == 0]
    area = [y * s for y, s in zip(ymax, sgma)]
    
    # Lines 11623-11625: Create and return result structure
    names = ['xMaxDeriv', 'xMaxGauss', 'sigma', 'yMaxFit', 'Area', 'ysm@xpeaks', 'aer-1*sigma', 'aer-2*sigma', 'aer+1*sigma', 'aer+2*sigma']
    PeakList = list(zip([x[p] for p in PEAKS], xpeaks, sgma, ymax, area, ymax2, aerM1s, aerM2s, aerP1s, aerP2s))
    return {'names': names, 'data': PeakList}


def PS2_ccc(x, y):
    # Line 11629: Set separation limit
    sepp = 8  # limit to separate independent peaks, unit: bins of bin mass scale, i.e. 64 ppm above 125 Da

    # Lines 11630-11631: Smooth y and calculate derivative
    ysm = smooth(y, 5)
    dysm = deriv(x, ysm)

    # Lines 11632-11633: Get length and create index array
    length = len(y)
    index = range(length)

    # Lines 11634-11638: Find peaks
    PEAKS = np.maximum.reduce([ysm, np.pad(ysm, (0, 2))[:length], np.pad(ysm, (1, 1))[:length], 
                               np.pad(ysm, (2, 0))[:length]])
    PEAKS[dysm < 0] = 0
    PEAKS[np.pad(dysm, (1, 0))[:length] > 0] = 0
    PEAKS[PEAKS < np.max(ysm) * getpar('LowThresh')] = 0
    PEAKS[PEAKS < 0.55] = 0
    PEAKS = [i for i, peak in enumerate(PEAKS) if peak > 0]

    # Line 11639: Get number of peaks
    length2 = len(PEAKS)

    print('HHGHGH', length2)

    # Lines 11641-11650: Initialize variables
    xpeaks, sgma, ymax, ymax2, aerM1s, aerM2s, aerP1s, aerP2s = [0], [0], [0], [0], [0], [0], [0], [0]
    notres = [1]

    # Lines 11651-11714: Main peak processing loop
    for i in range(length2):
        indis = range(PEAKS[i] - 5, PEAKS[i] + 6)
        xi = [x[j] for j in indis]
        yi = [y[j] for j in indis]

        # Gaussian fitting
        if max(yi) <= 10:
            coeff = gaussfit(xi, yi, nterms=3)
        else:
            coeff = gaussfit(xi, [y + 0.0001 for y in yi])

        if coeff[0] < 0 or coeff[0] > max(yi) * 3 or 1e6 * abs(coeff[2]) / np.mean(xi) < 5 or coeff[2] <= 0:
            coeff = gaussfit(xi, [ysm[j] for j in indis], nterms=3)

        # Calculate separation limit and other parameters
        seplim = max(float(sepp) * 1e6 * (x[np.argmin(np.abs(np.array(x) - coeff[1])) + 1] - 
                                          x[np.argmin(np.abs(np.array(x) - coeff[1]))]) / coeff[1])
        xpeaks.append(coeff[1])
        ymax.append(coeff[0])
        flt = np.argmin(np.abs(np.array(x) - coeff[1]))
        ymax2.append(ysm[flt])
        sgma.append(coeff[2])

        # Calculate aerM1s, aerM2s, aerP1s, aerP2s
        flt2 = [i for i, xi in enumerate(x) if coeff[1] - 2 * coeff[2] < xi <= coeff[1] - coeff[2]]
        aerM2s.append(sum([y[i] for i in flt2]))
        flt2 = [i for i, xi in enumerate(x) if coeff[1] - coeff[2] < xi <= coeff[1]]
        aerM1s.append(sum([y[i] for i in flt2]))
        flt2 = [i for i, xi in enumerate(x) if coeff[1] < xi <= coeff[1] + coeff[2]]
        aerP1s.append(sum([y[i] for i in flt2]))
        flt2 = [i for i, xi in enumerate(x) if coeff[1] + coeff[2] < xi <= coeff[1] + 2 * coeff[2]]
        aerP2s.append(sum([y[i] for i in flt2]))

        # Determine if peak is resolved
        if coeff[0] < 0:
            notres.append(1)
        elif coeff[0] > 2.5 * max(yi):
            notres.append(1)
        elif coeff[1] < min(xi) or coeff[1] > max(xi):
            notres.append(1)
        elif coeff[2] < 0:
            notres.append(1)
        elif min(abs(np.array(xpeaks[:-1]) - coeff[1]) / coeff[1]) < seplim * 1e-6:
            check = 1
            iii = [j for j, xp in enumerate(xpeaks[:-1]) if 1e6 * abs((xp - coeff[1]) / coeff[1]) < seplim]
            if min([notres[j] for j in iii]) == 1:
                notres.append(0)
            else:
                iv = [j for j in iii if notres[j] == 0]
                if max([ymax2[j] for j in iv]) > ysm[flt]:
                    notres.append(1)
                else:
                    for j in iv:
                        notres[j] = 1
                    notres.append(0)
        else:
            notres.append(0)

    # Lines 11716-11722: Filter results and create output structure
    valid_indices = [i for i, nr in enumerate(notres) if nr == 0]
    xpeaks = [xpeaks[i] for i in valid_indices]
    sgma = [sgma[i] for i in valid_indices]
    ymax = [ymax[i] for i in valid_indices]
    ymax2 = [ymax2[i] for i in valid_indices]
    aerM1s = [aerM1s[i] for i in valid_indices]
    aerM2s = [aerM2s[i] for i in valid_indices]
    aerP1s = [aerP1s[i] for i in valid_indices]
    aerP2s = [aerP2s[i] for i in valid_indices]
    PEAKS = [0] + [PEAKS[i-1] for i in valid_indices]
    area = [y * s for y, s in zip(ymax, sgma)]

    names = ['xMaxDeriv', 'xMaxGauss', 'sigma', 'yMaxFit', 'Area', 'ysm@xpeaks', 'aer-1*sigma', 'aer-2*sigma', 'aer+1*sigma', 'aer+2*sigma']
    PeakList = list(zip(PEAKS, xpeaks, sgma, ymax, area, ymax2, aerM1s, aerM2s, aerP1s, aerP2s))

    return {'names': names, 'data': PeakList}


def SmoothSumSpec_old(SumSpectrum, a, t0, ex, SampInt, instrument):
    # Line 11726: Set compilation options (not needed in Python)
    
    # Line 11730: Calculate factor based on SampInt
    faktor = max(1, int(1e-10 / SampInt))
    
    # Lines 11732-11733: Set f2 based on instrument type
    f2 = 2 if instrument.startswith('TOF1000') else 1
    
    # Lines 11735-11744: Apply Savitzky-Golay filters
    savgolFilter = savgol(9*faktor, 9*faktor, 0, int(10/(f2*faktor)))
    dataX0 = convolve(SumSpectrum, savgolFilter, mode='same')
    derdataX0 = convolve(np.gradient(dataX0), savgolFilter, mode='same')
    
    savgolFilter = savgol(9*faktor, 9*faktor, 0, int(6/(f2*faktor)))
    dataX1 = convolve(SumSpectrum, savgolFilter, mode='same')
    derdataX1 = convolve(np.gradient(dataX1), savgolFilter, mode='same')
    
    savgolFilter = savgol(14*faktor, 14*faktor, 0, int(6/(f2*faktor)))
    dataX2 = convolve(SumSpectrum, savgolFilter, mode='same')
    derdataX2 = convolve(np.gradient(dataX2), savgolFilter, mode='same')
    
    savgolFilter = savgol(20*faktor, 20*faktor, 0, int(6/(f2*faktor)))
    dataX3 = convolve(SumSpectrum, savgolFilter, mode='same')
    derdataX3 = convolve(np.gradient(dataX3), savgolFilter, mode='same')
    
    # Line 11745: Get number of points
    Npts = len(SumSpectrum)
    
    # Lines 11746-11748: Calculate thresholds
    Thresh1 = m2t(64.6, a, t0, ex, SampInt)
    Thresh2 = m2t(144.6, a, t0, ex, SampInt)
    Thresh3 = m2t(324.6, a, t0, ex, SampInt)
    
    # Lines 11749-11758: Combine smoothed data
    data2 = SumSpectrum.copy()
    if Npts > Thresh1:
        data2[:Thresh1] = dataX0[:Thresh1]
    else:
        data2 = dataX0
    
    if Npts > Thresh2:
        data2[Thresh1:Thresh2] = dataX1[Thresh1:Thresh2]
    elif Npts > Thresh1:
        data2[Thresh1:] = dataX1[Thresh1:]
    
    if Npts > Thresh3:
        data2[Thresh2:Thresh3] = dataX2[Thresh2:Thresh3]
        data2[Thresh3:] = dataX3[Thresh3:]
    elif Npts > Thresh2:
        data2[Thresh2:] = dataX2[Thresh2:]
    
    # Line 11759: Calculate derivative
    derivData = np.gradient(data2)
    
    # Lines 11761-11770: Combine smoothed derivatives
    if Npts > Thresh1:
        derivData[:Thresh1] = derdataX0[:Thresh1]
    else:
        derivData = derdataX0
    
    if Npts > Thresh2:
        derivData[Thresh1:Thresh2] = derdataX1[Thresh1:Thresh2]
    elif Npts > Thresh1:
        derivData[Thresh1:] = derdataX1[Thresh1:]
    
    if Npts > Thresh3:
        derivData[Thresh2:Thresh3] = derdataX2[Thresh2:Thresh3]
        derivData[Thresh3:] = derdataX3[Thresh3:]
    elif Npts > Thresh2:
        derivData[Thresh2:] = derdataX2[Thresh2:]
    
    # Line 11772: Return results as a dictionary
    return {'Data': data2, 'derivData': derivData, 'Npts': Npts}

def SmoothSumSpec(SumSpectrum, a, t0, ex, SampInt, instrument):
    # Line 11775: Get dimensions of SumSpectrum
    Npts = SumSpectrum.shape[0]
    
    # Lines 11776-11777: Set resolution based on instrument
    if instrument.startswith('TOF1000'):
        res = 2500
    elif instrument.startswith('VOCUS'):
        res = 10000
    else:
        res = 5000
    
    # Line 11783: Get smoothing factor from parameters
    factor = getpar('SmFact')
    
    # Lines 11784-11791: Calculate extension values for different mass ranges
    mm_values = [20.6, 64.6, 144.6, 324.6]
    ext_values = []
    for mm in mm_values:
        ext = floor(factor * (m2t(mm + mm/res, a, t0, 0.5, SampInt) - m2t(mm, a, t0, 0.5, SampInt))) + 1
        ext_values.append(ext)
    ext0, ext1, ext2, ext3 = ext_values
    
    # Lines 11793-11810: Apply Savitzky-Golay filters and convolutions
    for ext in ext_values:
        savgolFilter = savgol(ext, ext, 0, min([2*ext, 4+max(floor(factor/1.5))]))
        dataX = convolve(SumSpectrum, savgolFilter, mode='same')
        savgolFilter = savgol(ext, ext, 1, min([2*ext, 4+max(floor(factor/1.5))]))
        derdataX = convolve(SumSpectrum, savgolFilter, mode='same')
        globals()[f'dataX{ext_values.index(ext)}'] = dataX
        globals()[f'derdataX{ext_values.index(ext)}'] = derdataX
    
    # Lines 11814-11817: Calculate threshold values
    Thresh1 = m2t(64.6, a, t0, ex, SampInt)
    Thresh2 = m2t(144.6, a, t0, ex, SampInt)
    Thresh3 = m2t(324.6, a, t0, ex, SampInt)
    data2 = SumSpectrum.copy()
    
    # Lines 11818-11826: Apply smoothed data based on thresholds
    if Npts > Thresh1:
        data2[:Thresh1] = dataX0[:Thresh1]
    else:
        data2 = dataX0
    
    if Npts > Thresh2:
        data2[Thresh1:Thresh2] = dataX1[Thresh1:Thresh2]
    elif Npts > Thresh1:
        data2[Thresh1:] = dataX1[Thresh1:]
    
    if Npts > Thresh3:
        data2[Thresh2:Thresh3] = dataX2[Thresh2:Thresh3]
        data2[Thresh3:] = dataX3[Thresh3:]
    elif Npts > Thresh2:
        data2[Thresh2:] = dataX2[Thresh2:]
    
    # Line 11827: Calculate derivative of smoothed data
    derivData = np.gradient(data2)
    
    # Lines 11829-11838: Apply smoothed derivatives based on thresholds
    if Npts > Thresh1:
        derivData[:Thresh1] = derdataX0[:Thresh1]
    else:
        derivData = derdataX0
    
    if Npts > Thresh2:
        derivData[Thresh1:Thresh2] = derdataX1[Thresh1:Thresh2]
    elif Npts > Thresh1:
        derivData[Thresh1:] = derdataX1[Thresh1:]
    
    if Npts > Thresh3:
        derivData[Thresh2:Thresh3] = derdataX2[Thresh2:Thresh3]
        derivData[Thresh3:] = derdataX3[Thresh3:]
    elif Npts > Thresh2:
        derivData[Thresh2:] = derdataX2[Thresh2:]
    
    # Lines 11841: Return results as a dictionary
    return {'Data': data2, 'derivData': derivData, 'Npts': Npts}


def PeakTable(masslist: List[float], resolution: float = 4000) -> List[List[float]]:
    # Line 25365: Check if resolution is provided, otherwise use default
    if resolution is None:
        resolution = 4000

    # Line 25366-25367: Initialize IntList
    length: int = len(masslist)
    IntList: List[List[float]] = [[0.0] * 6 for _ in range(length)]

    # Line 25368: Set first column of IntList to masslist
    for i in range(length):
        IntList[i][0] = masslist[i]

    # Line 25369: Set resolution
    res: float = resolution

    # Lines 25370-25398: Main loop to calculate peak boundaries
    for i in range(length):
        # Line 25371: Calculate C (constant for peak width)
        C: float = 2 * (2 * math.log(2)) ** 0.5

        # Line 25372: Calculate sigma (peak width)
        sig: float = IntList[i][0] / (res * C)

        # Lines 25373-25375: Calculate peak boundaries and sigma
        IntList[i][1] = IntList[i][0] - 2 * sig
        IntList[i][2] = IntList[i][0] + 2 * sig
        IntList[i][3] = sig

        # Lines 25376-25377: Calculate Down and Up (adjacent peak boundaries)
        Down: float = IntList[i-1][0] * (1 + 2/(res*C)) if i > 0 else 0
        Up: float = IntList[i+1][0] * (1 - 2/(res*C)) if i < length - 1 else 100000

        # Lines 25378-25386: Adjust lower boundary if overlap with previous peak
        if Down > IntList[i][1]:
            overlap: float = Down - IntList[i][1]
            IntList[i][1] = IntList[i][0] - 2 * sig + overlap / 2
            IntList[i][4] = overlap / (2 * sig)

        # Lines 25387-25395: Adjust upper boundary if overlap with next peak
        if Up < IntList[i][2]:
            overlap: float = IntList[i][2] - Up
            IntList[i][2] = IntList[i][0] + 2 * sig - overlap / 2
            IntList[i][5] = overlap / (2 * sig)

    # Line 25399: Return the calculated IntList
    return IntList


def multiples(vec: List[float]) -> List[Tuple[float, int]]:
    """
    Find consecutive repeated values in a vector and count their occurrences.

    Args:
        vec (List[float]): Input vector of float values.

    Returns:
        List[Tuple[float, int]]: List of tuples containing unique values and their counts.
    """
    # Line 25343: Initialize result list
    res: List[Tuple[float, int]] = []
    
    # Line 25344: Get length of input vector
    length: int = len(vec)

    # Line 25347: Initialize index
    ii: int = 0

    # Lines 25349-25359: Main loop to find consecutive repeated values
    while ii < length:
        counter: int = 0
        # Check for consecutive repeated values
        while ii + counter + 1 < length and vec[ii + counter + 1] == vec[ii + counter]:
            counter += 1
            if ii + counter >= length - 1:
                break
        
        # Add the value and its count to the result
        res.append((vec[ii + counter], counter + 1))
        ii += counter + 1

    # Line 25361: Return the result
    return res


def ppm_bin(masslist: List[float], BinWidth_ppm: float) -> List[Tuple[float, float, float]]:
    # Line 25404: Convert BinWidth_ppm to fractional value
    BinWidth_ppm = BinWidth_ppm / 1e6

    # Lines 25406-25407: Set start values
    start = 10.0
    startppm = min(2000, 0.001 * int(1 / BinWidth_ppm))

    # Lines 25409-25410: Split masslist into ppm and mDa ranges
    m_ppm = [m for m in masslist if m >= startppm]
    m_mDa = [m for m in masslist if start < m < startppm]

    # Lines 25414-25425: Calculate ppm bins
    if m_ppm:
        m_start = [startppm * math.exp(BinWidth_ppm * i) for i in range(2000000)]
        m_end = m_start[1:] + [max(m_start) + 1]
        filt = [i for i, m in enumerate(m_start) if m < 2000]
        m_start2 = [m_start[i] for i in filt]
        m_end2 = [m_end[i] for i in filt]

        length = len(m_start2)
        scores2 = [0.0] * length
        I = [int(math.log(m / startppm) / BinWidth_ppm) for m in m_ppm]
        I.sort()
        I2 = multiples(I)
        for idx, count in I2:
            scores2[idx] = count

    # Lines 25428-25438: Calculate mDa bins
    if m_mDa:
        mDa_bins = int(1000 * (startppm - start))
        m_start = [start + i / 1000 for i in range(mDa_bins)]
        m_end = m_start[1:] + [max(m_start) + 0.001]

        length = len(m_start)
        scores = [0.0] * length
        I = [int((m - start) / 0.001) for m in m_mDa]
        I.sort()
        I2 = multiples(I)
        for idx, count in I2:
            scores[idx] = count

    # Lines 25439-25447: Combine results
    if not m_mDa:
        m_start, m_end, scores = m_start2, m_end2, scores2
    elif m_mDa and m_ppm:
        m_start = m_start + m_start2
        m_end = m_end + m_end2
        scores = scores + scores2

    # Line 25449: Return results
    return list(zip(m_start, m_end, scores))


def quantile_int(x: List[float], p: float, type: int) -> float:
    """
    Calculate sample quantiles.
    
    Args:
    x: List of float values to compute quantiles from.
    p: Probability for which to compute the quantile (0 <= p <= 1).
    type: Integer specifying the quantile algorithm (1-9).
    
    Returns:
    The computed quantile value.
    """
    # Lines 26307-26316: Machine precision setup
    eps = 2.2204460492503131e-16  # Double precision machine epsilon
    
    n_x = len(x)
    
    # Lines 26318-26348: Discrete estimators (types 1-3)
    if type <= 3:
        m_types = [0, 0, -0.5]
        m = m_types[type - 1]
        
        j = int(p * n_x + m)
        g = p * n_x + m - j
        
        if type == 1:
            gamma = 1.0 if g > 0 else 0.0
        elif type == 2:
            gamma = 1.0 if g > 0 else 0.5
        else:  # type == 3
            if g > 0:
                gamma = 1.0
            else:
                gamma = 0.0 if j % 2 == 0 else 1.0
    
    # Lines 26350-26365: Continuous estimators (types 4-9)
    else:
        alpha_all = [0, 0.5, 0., 1, 0.333333333, 0.375]
        alpha = alpha_all[type - 4]
        beta = alpha
        
        m = alpha + p * (1.0 - alpha - beta)
        j = int(p * n_x + m + 4.0 * eps)
        g = p * n_x + m - j
        gamma = 0.0 if abs(g) < 4 * eps else g
    
    # Lines 26367-26371: Calculate the quantile
    j = j - 1
    Q_p = (1.0 - gamma) * x[max(0, j)] + gamma * x[min(j + 1, n_x - 1)]
    
    return Q_p



def quantile(x: List[float], p: Union[float, List[float]], type: int = 8, sorted: bool = False, nan: bool = False) -> Union[float, List[float]]:
    """
    Calculate the quantiles of a distribution of points.

    Args:
        x (List[float]): The data vector.
        p (Union[float, List[float]]): The quantile levels requested (between 0 and 1).
        type (int, optional): The type of quantile calculation. Defaults to 8.
        sorted (bool, optional): If the data is already sorted. Defaults to False.
        nan (bool, optional): Filter out NaN values. Defaults to False.

    Returns:
        Union[float, List[float]]: The calculated quantile(s).
    """
    # Lines 26392-26395: Check type validity
    if type < 1 or type > 9:
        print('QUANTILE: Only type 1-9 is supported!')
        return -1

    # Lines 26401-26403: Determine if p is scalar or vector
    n_p = len(p) if isinstance(p, list) else 1
    is_scalar = isinstance(p, float)
    res = [0.0] * n_p

    # Lines 26405-26408: Sort data if necessary
    xs = x if sorted else sorted(x)

    # Lines 26410-26422: Handle NaN values if specified
    if nan:
        xs = [val for val in xs if not math.isnan(val)]
        if not xs:
            return float('nan') if is_scalar else [float('nan')] * n_p

    # Lines 26426-26427: Calculate quantiles
    for i in range(n_p):
        res[i] = quantile_int(xs, p[i] if isinstance(p, list) else p, type)

    # Lines 26429-26430: Return result
    return res[0] if is_scalar else res

# Helper functions (to be implemented)
def smooth_sum_spec(sum_spectrum, a, t0, ex, samp_int, instrument):
    # Implement smoothing logic here
    pass

def find_peak_characteristics(deriv_data, i, step, thresh):
    # Implement peak characteristic finding logic here
    pass

def clean_peak(data2, data3, i_peak_max, m_unit, peak_broadness, counts_max):
    # Implement peak cleaning logic here
    pass


# Helper functions (to be implemented)
def m2t(peaklist, a, t0, ex, sampint, time=False):
    # Implement m2t function
    pass



def optimize_exponent(M1, M2, M3, Const):
    # Implement exponent optimization
    pass


def generate_mass_variations(M1, M2, M3, var_ppm):
    # Implement mass variation generation
    pass
