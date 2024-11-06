from typing import List, Optional, Tuple, Union

def fitTR(m19: float, m37: float, MtoQvalues: List[float], expPPBs: List[float], 
          Sensitivities: List[float], Counts: List[float], BGcounts: Optional[List[float]] = None, 
          Figure: bool = False) -> List[float]:
    # Line 14882-14883: Initialize BGcounts if not provided
    if BGcounts is None:
        BGcounts = [0] * len(Counts)

    # Line 14884-14887: Initialize variables
    numMtoQ = len(MtoQvalues)
    indis = list(range(numMtoQ))
    NC = [0.0] * numMtoQ
    MeasSens = [0.0] * numMtoQ
    Set1 = [0.0] * numMtoQ

    # Line 14890-14899: Set initial parameters
    mLOW, wLOW = 10, 2
    mHIGH, wHIGH = 1000, 10
    ex = 0.5
    DELTAex = 0.01
    EXrange = [0, 2]
    field = [0.0] * 201
    scoreee = 1e10
    bestEX = 0.5

    # Line 14922-14970: Main loop to find best EX
    for dEX in range(201):
        ex = DELTAex * dEX
        
        # Calculate NC and MeasSens
        for i in range(numMtoQ):
            factor1 = tr(MtoQvalues[i], ex, mLOW, wLOW, mHIGH, wHIGH) / tr(19, ex, mLOW, wLOW, mHIGH, wHIGH)
            factor2 = tr(MtoQvalues[i], ex, mLOW, wLOW, mHIGH, wHIGH) / tr(37, ex, mLOW, wLOW, mHIGH, wHIGH)
            NC[i] = (Counts[i] - BGcounts[i]) * 1e6 / (factor1 * m19 + factor2 * m37)

        for i in range(numMtoQ):
            MeasSens[i] = NC[i] / expPPBs[i]
            # Add contributions from other masses
            if abs(MtoQvalues[i] - 69) < 0.5:
                MeasSens[i] += NC[next(j for j, m in enumerate(MtoQvalues) if abs(m - 41) < 0.5)] / expPPBs[i]
            # ... (similar adjustments for other masses)
            Set1[i] = MeasSens[i] / Sensitivities[i]

        # Filter and process data
        fiiil = [i for i, m in enumerate(MtoQvalues) if 58.5 < m < 138.5 and Set1[i] < 9999]
        set_values = [Set1[i] for i in fiiil]
        xx = [MtoQvalues[i] for i in fiiil]
        set_values = [s / max(set_values) for s in set_values]

        # Calculate slope and update bestEX
        slpe = abs(max(s for s, x in zip(set_values, xx) if x < 74) - 
                   max(s for s, x in zip(set_values, xx) if x > 120))
        if slpe < scoreee:
            scoreee = slpe
            bestEX = ex

        # Plotting code would go here (omitted for brevity)

    # Line 14972-15082: Iterative optimization of bestEX
    bestEXold = -9999
    while abs(bestEXold - bestEX) > 0.02:
        # ... (optimization logic)
        pass  # Placeholder for the optimization loop

    # Line 15084-15124: Find best HIGH filter
    # ... (HIGH filter optimization logic)

    TRpar = [bestex, mlow, wlow, mhigh, whigh]

    # Line 15131-15167: Generate figure if requested
    if Figure:
        # ... (figure generation logic)
        pass  # Placeholder for figure generation

    return TRpar

def caNC(ex: float, mLOW: float, wLOW: float, mHIGH: float, wHIGH: float, 
         m19: float, m37: float, MtoQvalues: List[float], Counts: List[float], 
         BGcounts: Union[List[float], None] = None, 
         low: bool = False, medium: bool = False, high: bool = False) -> List[float]:
    """
    Calculate normalized counts (NC) based on input parameters.
    
    # Line 15170: Function definition with keyword arguments
    """
    
    # Line 15172: Set BGcounts to zero if not provided
    if BGcounts is None:
        BGcounts = [0] * len(Counts)
    
    # Lines 15174-15176: Initialize variables
    numMtoQ = len(MtoQvalues)
    indis = list(range(numMtoQ))
    NC = [0.0] * numMtoQ
    
    # Lines 15178-15186: Filter indices based on keyword arguments
    if low:
        indis = [i for i, v in enumerate(MtoQvalues) if v < 80.5]
    if medium:
        indis = [i for i, v in enumerate(MtoQvalues) if 58.5 < v < 138.5]
    if high:
        indis = [i for i, v in enumerate(MtoQvalues) if v > 106.5]
    
    # Lines 15188-15193: Calculate NC for each index
    for i in indis:
        factor1 = tr(MtoQvalues[i], ex, mLOW, wLOW, mHIGH, wHIGH) / tr(19, ex, mLOW, wLOW, mHIGH, wHIGH)
        factor2 = tr(MtoQvalues[i], ex, mLOW, wLOW, mHIGH, wHIGH) / tr(37, ex, mLOW, wLOW, mHIGH, wHIGH)
        NC[i] = (Counts[i] - BGcounts[i]) * 1e6 / (factor1 * m19 + factor2 * m37)
    
    # Line 15194: Return the calculated NC
    return NC


def calcNC(ex: float, masslow: float, widthlow: float, MASShigh: float, WIDTHhigh: float, 
           counts: List[List[float]], countsBG: List[List[float]], ADD40: int, 
           low: Optional[bool] = False, medium: Optional[bool] = False, high: Optional[bool] = False) -> List[float]:
    # Line 15197-15206: Initialize indices based on keyword arguments
    indis = list(range(22))
    if low:
        indis = [0, 1, 2, 3, 7, 8, 9]
    if medium:
        indis = [3, 5, 6, 7, 8, 9, 10, 11, 12, 13, 14]
    if high:
        indis = [10, 11, 12, 13, 14, 18, 19, 20, 21]

    # Line 15209: Define mass-to-charge ratios
    mTOz = [33, 42, 45, 59, 87, 69, 41, 71, 73, 79, 107, 121, 137, 81, 133, 181, 223, 207, 297, 281, 371, 355]

    # Lines 15210-15214: Initialize arrays
    factor1 = [0.0] * 22
    factor2 = [0.0] * 22
    primi = [0.0] * 22
    C = [0.0] * 22
    NC = [0.0] * 22

    # Line 15217-15230: Main calculation loop
    for i in range(len(indis)):
        idx = indis[i]
        # Calculate transmission factors
        factor1[idx] = tr(mTOz[idx], ex, masslow, widthlow, MASShigh, WIDTHhigh) / tr(19, ex, masslow, widthlow, MASShigh, WIDTHhigh)
        factor2[idx] = tr(mTOz[idx], ex, masslow, widthlow, MASShigh, WIDTHhigh) / tr(37, ex, masslow, widthlow, MASShigh, WIDTHhigh)
        
        # Calculate primary ion signal
        primi[idx] = (sum(countsBG[ADD40:ADD40+10][28]) * factor1[idx] + sum(countsBG[ADD40:ADD40+10][29]) * factor2[idx]) / 4
        
        # Special cases for certain indices
        if idx == 14 or idx == 15 or idx == 9:
            primi[idx] = (sum(countsBG[ADD40:ADD40+10][28]) * factor1[idx] + 1 * sum(countsBG[ADD40:ADD40+10][29]) * factor2[idx]) / 4
        
        # Calculate total signal and normalized counts
        C[idx] = sum(counts[ADD40:ADD40+10][idx]) - sum(countsBG[ADD40:ADD40+10][idx])
        NC[idx] = C[idx] * 1e6 / primi[idx]

    return NC

# Helper function (needs to be implemented based on IDL's tr function)
def tr(m: float, ex: float, masslow: float, widthlow: float, MASShigh: float, WIDTHhigh: float) -> float:
    # Implement the transmission function logic here
    pass


def tr( float, ex: float, mLOW: float, wLOW: float, mHIGH: float, wHIGH: float, norm: bool = False) -> float:
    # Transmission function implementation
    # This function needs to be implemented based on the specific requirements
    pass

def smooth(data: List[float], window_size: int, edge_mirror: bool = False) -> List[float]:
    # Smoothing function implementation
    # This function needs to be implemented based on the specific requirements
    pass



def testSCALE(masses: List[float], lib: List[float], instrument: str, ions: List[float], devs: List[float]) -> dict[str, Union[float, List[float]]]:
    # Initialize parameters
    ppm: float = 20.0
    scor: float = 0.0
    scorppm: float = 0.0
    accuracy: float = 0.0
    LowDa: float = 0.0015

    # Adjust parameters based on instrument type (Lines 25831-25836)
    if instrument.startswith('TOF1000'):
        ppm = 100.0
        LowDa = 0.0075
    elif instrument.startswith('TOF8000') or instrument.startswith('VOCUS'):
        ppm = 20.0
        LowDa = 0.0015

    # Calculate initial deviation (Lines 25844-25845)
    length: int = len(masses)
    deviation: float = lib[min(range(len(lib)), key=lambda i: abs(lib[i] - masses[0]))] - masses[0]

    # Check for low mass matches (Lines 25847-25853)
    MatchLow: int = sum(1 for m in [21.0221, 59.0491, 30.994, 33.9935, 47.9966] if min(abs(masses - m)) < LowDa)
    if MatchLow < 1.5:
        deviation = 9999

    # Calculate deviations for different mass ranges (Lines 25861-25873)
    if deviation < 9999:
        split: List[float] = [masses[0], 120.6, 250.6, max(masses)]
        for sss in range(3):
            massi: List[float] = [m for m in masses if split[sss] < m <= split[sss+1]]
            libi: List[float] = [l for l in lib if split[sss] < l <= split[sss+1]]
            for m in massi:
                deviation = [deviation] + [libi[min(range(len(libi)), key=lambda i: abs(libi[i] - m))] - m]
            if sss == 0:
                devi120 = deviation

    # Calculate scores (Lines 25885-25893)
    filter: List[int] = [i for i, d in enumerate(deviation) if abs(d) < 0.002]
    filter2: List[int] = [i for i, d in enumerate(deviation) if abs(d) < 0.010]
    filterppm: List[int] = [i for i, d in enumerate(deviation) if 1e6 * abs(d) / masses[i] < ppm]

    if filter:
        scor = float(len(filter)) + float(len(filter2)) / 1000
    if filterppm:
        scorppm = float(len(filterppm)) + 0.2 - sum(1e6 * abs(deviation[i]) / masses[i] for i in filterppm) / (100 * len(filterppm))

    if filter:
        accuracy = sum(deviation[i] for i in filter) / scor

    # Adjust scores based on specific ions (Lines 25895-25900)
    for i in range(3):
        if ions[i] > 0 and min(abs(masses - ions[i])) > devs[i]:
            scor /= 3
            scorppm /= 3

    # Create and return result structure (Lines 25902-25904)
    s1: dict[str, Union[float, List[float]]] = {
        'scor': scor,
        'scorppm': scorppm,
        'accuracy': accuracy,
        'deviation': deviation,
        'masses': masses
    }
    return s1



def testSCALE2(massesAll: List[float], lib: List[float], instrument: str, ions: List[float], devs: List[float]) -> dict[str, Union[float, List[float]]]:
    # Lines 25910-25918: Set instrument-specific parameters
    ppm = 20.001
    LowDa = 0.0015
    if 'TOF1000' in instrument:
        ppm, LowDa = 100.0, 0.0075
    elif 'TOF8000' in instrument or 'VOCUS' in instrument:
        ppm, LowDa = 20.0, 0.0015

    # Lines 25921-25925: Initialize scoring variables
    scor, scor120, rat120 = 0.00001, 0.00001, 1.00001
    massesAll = [m for m in massesAll if m > 20]

    # Lines 25928-25932: Filter masses based on mass defect
    massesAll = [m for m in massesAll if 0.9*m/800 - 0.25 < (m - int(m + 0.3)) < 0.9*m/800 + 0.025]

    # Lines 25937-25947: Handle H3O+ mode specific filtering
    if all(int(l + 0.4) % 2 != 0 for l in lib if l > 50):
        m34 = [m for m in massesAll if 33.9 < m < 34.1]
        m48 = [m for m in massesAll if 47.9 < m < 48.1]
        massesAll = [m for m in massesAll if int(m + 0.4) % 2 != 0]
        massesAll.extend(m34 + m48)

    # Lines 25952-25967: Further filtering for large mass lists
    if len(massesAll) > 3000:
        masses = [m for m in massesAll if any(low < m < high for low, high in 
                  [(50, 65), (100, 115), (150, 165), (200, 215), (250, 265), (300, 315), (350, 365), (400, 415)])]
    else:
        masses = massesAll

    # Lines 25971-25981: Initial scoring and low mass matching
    deviation = min(abs(l - masses[0]) for l in lib)
    if deviation < ppm * 1e-6 * masses[0]:
        scor += 1
    MatchLow = sum(1 for m in [21.0221, 59.0491, 30.994, 33.9935, 47.9966] if min(abs(ma - m) for ma in massesAll) < LowDa)
    if MatchLow < 1.5:
        deviation = 9999

    # Lines 25988-26013: Main scoring loop
    if deviation < 9999:
        split = [20.6, 120.6, 250.6, max(masses)]
        for sss in range(3):
            massi = [m for m in masses if split[sss] < m <= split[sss+1]]
            libi = [l for l in lib if split[sss] < l <= split[sss+1]]
            for m in massi:
                deviation = min(abs(l - m) for l in libi)
                if deviation < ppm * 1e-6 * m or deviation < 0.0015:
                    scor += 1
                    devi = 1e4 * deviation / m
            if sss == 0:
                scor120, rat120 = scor, scor / len(massi)

    # Lines 26016-26023: Final score adjustments
    if scor > 0.5:
        scor += 0.2 - devi / scor
    if scor120 > 0.5:
        scor120 += 0.2 - devi / scor120

    # Lines 26025-26030: Ion-specific score adjustments
    for ion, dev in zip(ions, devs):
        if ion > 0 and min(abs(m - ion) for m in massesAll) > dev:
            scor /= 3
            scor120 /= 3

    # Lines 26031-26034: Final score check
    if rat120 < 0.35:
        scor, scor120 = 0.001, 0.001

    # Lines 26040-26041: Return results
    return {'scor': scor, 'scor120': scor120, 'masses': masses, 'rat120': rat120}


def testSCALE3(masses: List[float], lib: List[float], instrument: str) -> float:
    # Find specific masses of interest
    m21 = masses[min(range(len(masses)), key=lambda i: abs(masses[i] - 21.022))]
    m30 = masses[min(range(len(masses)), key=lambda i: abs(masses[i] - 29.997))]
    m32 = masses[min(range(len(masses)), key=lambda i: abs(masses[i] - 31.989))]
    m39 = masses[min(range(len(masses)), key=lambda i: abs(masses[i] - 39.033))]

    # Combine specific masses with masses between 40 and 250
    masses = [m21, m30, m32, m39] + [m for m in masses if 40 < m < 250]
    length = len(masses)

    # Calculate deviations from library values
    deviation = [min(lib, key=lambda x: abs(x - m)) - m for m in masses]

    # Set ppm threshold based on instrument type
    ppm = 50
    if instrument.startswith('TOF1000'):
        ppm = 200
    elif instrument.startswith('TOF8000') or instrument.startswith('VOCUS'):
        ppm = 50

    # Filter masses within ppm threshold
    filterppm = [i for i, d in enumerate(deviation) if abs(d) / masses[i] * 1e6 < ppm]

    # Calculate score if filtered masses exist
    if filterppm:
        scorppm = float(len(filterppm)) + 0.2 - sum(abs(deviation[i]) / masses[i] * 1e6 for i in filterppm) / (100 * len(filterppm))
    else:
        scorppm = 0

    # Calculate additional statistics (not used in final score)
    stats = [sorted(abs(m - round(m)) for m in masses)[i] for i in [int(0.3 * length), int(0.5 * length), int(0.7 * length)]]

    # Return inverted score
    return 1000 - scorppm

def var_exists(variable: Union[int, float, str, list, dict]) -> int:
    return 1 if variable is not None else -1
