from typing import List, Tuple, Dict, Any, Union
import numpy as np
import matplotlib.pyplot as plt

def ActrisProtocol_simple(path: str, name: str, ProtData: np.ndarray, masses: np.ndarray, 
                          ProtEngData: np.ndarray, engnames: List[str], computerID: str) -> Dict[str, Any]:
    """
    Perform ACTRIS protocol analysis on PTR-MS data.
    
    Args:
        path (str): Path to data directory
        name (str): Name of the data file
        ProtData (np.ndarray): Protocol data array
        masses (np.ndarray): Mass values
        ProtEngData (np.ndarray): Engineering data array
        engnames (List[str]): Names of engineering parameters
        computerID (str): Identifier for the computer/instrument

    Returns:
        Dict[str, Any]: Results of the analysis
    """
    # Lines 15234-15249: Set default parameters
    T_FC: float = 31
    Fcal: float = 200
    Fha2o: List[float] = [0.0, 0.0, 0.0]

    if max(ProtEngData[0]) > 3770:
        Fcal = 240.0 * 298 / 273

    Fplus: float = 0
    if 'VOCUS PTR-TOF' in computerID:
        Fcal = 2000
    if len(computerID) == 0:
        Fcal = 2000
        Fplus = 5000 * 298 / 273  # Harvard

    Flow: float = float(Fcal) * 273 / 298 + Fplus + 9.6

    # Lines 15251-15254: Calculate dilution and loop temperature
    Dil: float = (9.6 * 273 / (273 + T_FC)) / Flow
    T_loop: float = 273 + 100

    # Lines 15258-15262: Create output directory
    Path2: str = path + 'STD_PROT/'
    os.makedirs(Path2, exist_ok=True)

    # Lines 15263-15267: Calculate relative time and extract PTR ID
    starttim: float = min(ProtEngData[:, 0])
    reltime: np.ndarray = (ProtEngData[:, 0] - starttim) * 24 * 60
    PTRid: str = os.path.basename(os.path.dirname(Path2))

    # Lines 15268-15272: Convert cps to counts
    timestep: float = max(reltime) * 60 / len(reltime)

    # Lines 15275-15309: Extract and process ion data
    I19: np.ndarray = extract_ion_data(ProtData, masses, 19.017)
    I21: np.ndarray = extract_ion_data(ProtData, masses, 21.022)
    # ... (similar extractions for other ions)

    # Lines 15312-15316: Calculate primary ions
    Prim: np.ndarray = I19 + I37

    # Lines 15319-15321: Calculate ion fractions
    F37: np.ndarray = I37 / Prim
    F32: np.ndarray = I34 * 242 / Prim

    # Lines 15323-15336: Extract engineering data
    Pdrift: np.ndarray = extract_eng_data(ProtEngData, engnames, 'p-Drift')
    H2Oscc: np.ndarray = extract_eng_data(ProtEngData, engnames, 'H2O scc')
    # ... (similar extractions for other engineering data)

    # Lines 15338-15346: Adjust data for specific instruments
    if 'VOCUS PTR-TOF' in computerID:
        Pdrift = extract_eng_data(ProtEngData, engnames, 'SSQ pressure monitor [mbar]')
        Tdrift = extract_eng_data(ProtEngData, engnames, 'TOF temp monitor')
        Udrift = extract_eng_data(ProtEngData, engnames, 'Spare 1 monitor [V]') - 25

    # Lines 15349-15424: Extract and process compound data
    D5_355: np.ndarray = extract_compound_data(ProtData, masses, 355.07, 1.67)
    D5_371: np.ndarray = extract_compound_data(ProtData, masses, 371.102, 1.69)
    # ... (similar extractions for other compounds)

    # Lines 15427-15461: Calculate total signal
    TotSIG: np.ndarray = calculate_total_signal([METH_33, CH3CN_42, Acetal_45, Ace_59, MBO_87, 
                                                 MBO_69, MBO_41, MVK_71, MEK_73, Ben_79,
                                                 Xyl_107, TMB_121, Pin_137, Pin_81, benzF3_133, 
                                                 benzCl3_181, D3_223, D3_207, D4_297, D4_281,
                                                 D5_371, D5_355])

    # Lines 15464-15503: Extract isotope data
    D5_358iso: np.ndarray = extract_isotope_data(ProtData, masses, 358.073, 0.0358)
    D5_374iso: np.ndarray = extract_isotope_data(ProtData, masses, 374.105, 0.0370)
    # ... (similar extractions for other isotopes)

    # Lines 15506-15512: Calculate total isotope signal
    TotSIGiso: np.ndarray = calculate_total_signal([METH_34iso, CH3CN_43iso, Acetal_46iso, Ace_60iso, 
                                                    MBO_88iso, MBO_70iso, MBO_42iso, MVK_72iso, 
                                                    MEK_74iso, Ben_80iso, Xyl_108iso, TMB_122iso, 
                                                    Pin_138iso, Pin_82iso, benzF3_134iso, benzCl3_182iso, 
                                                    D3_225iso, D3_209iso, D4_299iso, D4_283iso,
                                                    D5_374iso, D5_358iso])

    # Lines 15515-15521: Define mass-to-charge ratios and compound names
    mTOz: List[float] = [33, 42, 45, 59, 87, 69, 41, 71, 73, 79, 107, 121, 137, 81, 133, 181, 223, 207, 297, 281,
                         371, 355, -1, -1, -1, -1, -1, -1, 19, 37, 39, 34, 43, 46, 60, 88, 70, 42, 72, 74,
                         80, 108, 122, 138, 82, 134, 182, 225, 209, 299, 283, 374, 358, -1, -1, -1, 78, 32, 120, 106]

    compounds: List[str] = ['METH_33', 'CH3CN_42', 'Acetal_45', 'Ace_59', 'MBO_87', 'ISO/MBO_69', 'ISO/MBO_41', 'MVK_71', 'MEK_73', 'Ben_79',
                            'Xyl_107', 'TMB_121', 'MT_137', 'MT_81', 'benzF3_133', 'benzCl3_181', 'D3_223', 'D3_207', 'D4_297', 'D4_281',
                            'D5_371', 'D5_355', 'TotSIG', 'FCinlet', 'reltime', 'TotSIGiso', 'prim', 'Tdrift', 'I19', 'I37',
                            'I39', 'METH_34iso', 'CH3CN_43iso', 'Acetal_46iso', 'Ace_60iso', 'MBO_88iso', 'MBO_70iso', 'MBO_42iso', 'MVK_72iso', 'MEK_74iso', 'Ben_80iso',
                            'Xyl_108iso', 'TMB_122iso', 'Pin_138iso', 'Pin_82iso', 'benzF3_134iso', 'benzCl3_182iso', 'D3_225iso', 'D3_209iso', 'D4_299iso', 'D4_283iso',
                            'D5_374iso', 'D5_358iso', 'Udx', 'pdrift', 'Udrift', 'BEN_78', 'I32', 'TMB_120', 'XYL_106']

    # Lines 15524-15527: Define isotope correction factors and reaction rate constants
    isotopecor: List[float] = [1.01, 1.02, 1.02, 1.03, 1.06, 1.06, 1.03, 1.05, 1.05, 1.07, 1.09, 1.1, 1.12, 1.07, 1.07, 2.48, 1.37, 1.36, 1.52, 1.51, 1.69, 1.67, 1, 1, 1, 1, 1, 1]
    RRK: List[float] = [2.2/1.3, 2.6+0.5, 3.03, 3.25, 2.68/1.45, 2.68/1.45, 2.68/1.45, 3.53/1.3, 3.25, 1.97, 2.31, 2.4, 2.45/1.2, 2.45/1.2, 2.46, 2.88/1.2, 2.59/1.2, 2.59/1.2, 2.99, 2.99, 3.39, 3.39*1.0, 3]
    RRK = [x * 1e-9 for x in RRK]

    # Lines 15529-15532: Define volume mixing ratios
    VMRppmNPL: List[float] = [1.019, 1.020, 1.001, 0.983, 0.996, 0.996, 0.996, 0.961, 1.009, 1.025, 0.998, 1.001, 0.989, 0.989, 1.047, 0, 0, 0, 0.901, 0.901, 1.051, 1.051, 4, 0, 0, 0, 0, 0]
    VMRppmAR: List[float] = [1.011, 1.01, 0.6, 0.967, 0.998, 0.998, 0.998, 0.55, 1.017, 1.006, 0.983, 0.989, 0.983, 0.983, 1.2, 0.995, 0.992, 0.992, 0.995, 0.995, 0.995, 0.995, 4, 0, 0, 0, 0, 0]

    # Lines 15534-15564: Combine all compound data
    CompARR: np.ndarray = np.array([METH_33, CH3CN_42, ACETAL_45, ACE_59, MBO_87,
                                    MBO_69, MBO_41, MVK_71, MEK_73, BEN_79,
                                    XYL_107, TMB_121, Pin_137, Pin_81, benzF3_133,
                                    benzCl3_181, D3_223, D3_207, D4_297, D4_281,
                                    D5_371, D5_355, TotSIG, np.full_like(FCinlet, Flow), reltime,
                                    TotSIGiso, Prim, Tdrift, I19, I37,
                                    I39, METH_34iso, CH3CN_43iso, Acetal_46iso, Ace_60iso, MBO_88iso, MBO_70iso, MBO_42iso, MVK_72iso, MEK_74iso,
                                    Ben_80iso, Xyl_108iso, TMB_122iso, Pin_138iso, Pin_82iso, benzF3_134iso, benzCl3_182iso, D3_225iso, D3_209iso, D4_299iso, D4_283iso,
                                    D5_374iso, D5_358iso, Udx, Pdrift, Udrift, BEN_78, I32, TMB_120, XYL_106])

    # Lines 15566-15569: Define indices for specific data
    i_Tdrift, i_Udx, i_pdrift, i_Udrift = 27, 53, 54, 55

    # Lines 15577-15706: Find injections (complex logic, simplified here)
    injections = find_injections(ACE_59, reltime, computerID)

    # Lines 15709-15764: Transmission correction (simplified)
    CompARRtranse, PrimTR, TotSigTR, countsTR, countsBGTR = transmission_correction(CompARR, mTOz, bestEX, MASSlow, WIDTHlow, MASShigh, WIDTHhigh)

    # Lines 15772-15830: Calculate sensitivities and other

    helpi = 1
    try:
        # Lines 16795-16801: Read HDF5 file attributes
        file_id = h5py.File(name, 'r')
        n_attr = len(file_id.attrs)
        ComputerID = 'nix'
        for j in range(n_attr):
            if 'Computer ID' in file_id.attrs.keys()[j]:
                ComputerID = file_id.attrs[file_id.attrs.keys()[j]]
        file_id.close()
    except:
        helpi = 0
        ComputerID = 'QMS_LSCE'

    # Lines 16805-16813: Set initial parameters
    T_FC = 31
    Fcal = 200
    if max(protengdata[0]) > 3770:
        Fcal = 240.0 * 298 / 273  # Orleans May 2019
    if 3800 < max(protengdata[0]) < 4005:
        Fcal = 240.0 * 291 / 273  # Sonnblick Aug-Dec 2019
        T_FC = 18  # Sonnblick Aug-Dec 2019

    # Lines 16815-16822: Adjust Fcal and Fplus based on ComputerID
    Fplus = 0
    if '2197a787-6c47-4c78-88d1-28fb6598bcc2' in ComputerID:
        Fcal = 2000  # VOCUS PTR-TOF
    if len(ComputerID) == 0:
        Fcal = 2000
        Fplus = 5000 * 298 / 273  # Harvard
    Flow = float(Fcal) * 273 / 298 + Fplus + 9.6

    # Lines 16825-16829: Calculate Dil and T_loop
    Dil = (9.6 * 273 / (273 + T_FC)) / Flow
    T_loop = 273 + 100
    if 3800 < max(protengdata[0]) < 4005:
        T_loop = 273 + 110  # Sonnblick Aug-Dec 2019

    # Lines 16831-16832: Adjust masses for specific ComputerID
    # This section is commented out in the original IDL code

    # Lines 16834-16839: Create output directory
    Path2 = os.path.join(path, 'STD_PROT/')
    if not os.path.exists(Path2):
        os.makedirs(Path2)

    # Lines 16840-16844: Calculate time-related variables
    starttim = min(ProtEngData[:, 0])
    reltime = (ProtEngData[:, 0] - min(ProtEngData[:, 0])) * 24 * 60
    hilfi = os.path.basename(os.path.dirname(Path2))
    PTRid = hilfi.split(os.path.sep)[0]
    Datumzeit = t09str(ProtEngData[0])  # Assuming t09str is defined elsewhere

    # Lines 16847-16848: Convert cps to counts
    timestep = max(reltime) * 60 / max(reltime.shape)

    # Lines 16851-16861: Extract specific ion data
    I19 = ProtData[np.abs(masses - 19.017).argmin(), :]
    I21 = ProtData[np.abs(masses - 21.022).argmin(), :]
    I37 = ProtData[np.abs(masses - 37.028).argmin(), :]
    I38 = ProtData[np.abs(masses - 38.033).argmin(), :]
    I39 = ProtData[np.abs(masses - 39.033).argmin(), :]
    I34 = ProtData[np.abs(masses - 33.994).argmin(), :]
    I32 = ProtData[np.abs(masses - 31.989).argmin(), :]
    I32[I32 > 6000] = I34[I32 > 6000] * 242

    # Lines 16863-16866: Adjust I19 and I37 based on conditions
    if max(I19) > 5000:
        I19 = I21 * 487
    I37[I37 < I38 * 645] = I38[I37 < I38 * 645] * 645

    # Lines 16869-16873: Calculate Prim, F37, F32
    Prim = I19 + I37
    F37 = I37 / Prim
    F32 = I34 * 242 / Prim

    # Lines 16875-16885: Extract engineering data
    Pdrift = ProtEngData[:, np.where(np.array([x.startswith('p-Drift') for x in engnames]))[0][0]]
    H2Oscc = ProtEngData[:, np.where(np.array([x.startswith('H2O scc') or x.startswith('H2O[scc') for x in engnames]))[0][0]]
    Udrift = ProtEngData[:, np.where(np.array([x.startswith('Udrift') for x in engnames]))[0][0]]
    Udx = ProtEngData[:, np.where(np.array([x.startswith('Udx') for x in engnames]))[0][0]]
    FCinlet = ProtEngData[:, np.where(np.array([x.startswith('FCinlet') for x in engnames]))[0][0]]
    Tdrift = ProtEngData[:, np.where(np.array([x.startswith('T-Drift') for x in engnames]))[0][0]]
    if '64db5058-a18e-4b09-ab32-f7cecac8afa4' in ComputerID:
        Tdrift = ProtEngData[:, -1] + 60  # Juelich PTR-TOF
    PC = ProtEngData[:, np.where(np.array([x.startswith('PC mba') or x.startswith('PC[mba') for x in engnames]))[0][0]]

    # Lines 16887-16892: Adjust parameters for VOCUS PTR-TOF
    if '2197a787-6c47-4c78-88d1-28fb6598bcc2' in ComputerID:
        Pdrift = ProtEngData[:, np.where(np.array([x.startswith('SSQ pressure monitor [mbar]') for x in engnames]))[0][0]]
        Tdrift = ProtEngData[:, np.where(np.array([x.startswith('TOF temp monitor') for x in engnames]))[0][0]]
        Udrift = ProtEngData[:, np.where(np.array([x.startswith('Spare 1 monitor [V]') for x in engnames]))[0][0]] - 25

    # Lines 16893-16897: Adjust parameters for PTR3
    if len(ComputerID) == 0:
        Tdrift = ProtEngData[:, np.where(np.array([x.startswith('TPS body temp monitor') for x in engnames]))[0][0]]
        Pdrift = ProtEngData[:, np.where(np.array([x.startswith('REAC   sensor: APR/CMR') for x in engnames]))[0][0]]
    
    # Lines 16898-16902: Set up parameters
    RelTime: np.ndarray = (engdata[:, 0] - min(engdata[:, 0])) * 24 * 60
    NameMarker: str = 'DO1'
    MarkerValue: float = 1.0
    DurationProtocol: float = 8.0

    # Lines 16904-16910: Find marker indices
    MarkerInd: np.ndarray = np.where(engdata[:, np.where(np.array([x.startswith(NameMarker) for x in engdatanames]))[0][0]] == MarkerValue)[0]
    StartInd: List[int] = [MarkerInd[0]]
    EndInd: List[int] = []
    for i in range(1, len(MarkerInd)):
        if MarkerInd[i] - MarkerInd[i-1] > 1:
            StartInd.append(MarkerInd[i])
            EndInd.append(MarkerInd[i-1])
    EndInd.append(MarkerInd[-1])

    # Lines 16912-16920: Process each section
    S_val: List[float] = []
    Serr_val: List[float] = []
    for i in range(len(StartInd)):
        start_time: float = RelTime[StartInd[i]]
        end_time: float = RelTime[EndInd[i]]
        if end_time - start_time > DurationProtocol:
            end_time = start_time + DurationProtocol
        
        # Lines 16922-16930: Calculate mean and standard deviation
        mean_cps: np.ndarray = np.mean(cps[StartInd[i]:EndInd[i], :], axis=0)
        std_cps: np.ndarray = np.std(cps[StartInd[i]:EndInd[i], :], axis=0)
        S_val.extend(mean_cps)
        Serr_val.extend(std_cps)

    # Lines 16932-16940: Calculate Fptr
    H3O: float = S_val[np.where(abs(masses - 21.0221) < 0.01)[0][0]]
    H3O18: float = S_val[np.where(abs(masses - 23.0274) < 0.01)[0][0]]
    NO: float = S_val[np.where(abs(masses - 30.9950) < 0.01)[0][0]]
    O2: float = S_val[np.where(abs(masses - 33.9935) < 0.01)[0][0]]
    Fptr: float = (H3O + H3O18 + NO + O2) / 1e6

    # Lines 16942-16950: Calculate additional parameters
    meanPdrift: float = np.mean(Pdrift[StartInd[0]:EndInd[-1]])
    meanTdrift: float = np.mean(Tdrift[StartInd[0]:EndInd[-1]])
    meanUdrift: float = np.mean(engdata[StartInd[0]:EndInd[-1], np.where(np.array([x.startswith('Udrift') for x in engdatanames]))[0][0]])

    S_val2 = [0.0] * 23
    Serr_val2 = [0.0] * 23

    for i in range(len(tStart)):
        mask = [(tStart[i] <= t <= tEnd[i]) for t in RelTime]
        cps_i = [row[mask] for row in cps]
        
        for j, (mass, row) in enumerate(zip(masses, cps_i)):
            if mass in [21.022, 39.033]:
                mean = sum(row) / len(row)
                std = ((sum((x - mean) ** 2 for x in row)) / (len(row) - 1)) ** 0.5
                S_val2[i] += mean
                Serr_val2[i] += std ** 2

    Fha2o[0] = S_val2[1] / S_val2[0] if S_val2[0] != 0 else 0
    Fha2o[1] = S_val[1] / S_val[0] if S_val[0] != 0 else 0
    Fha2o[2] = np.sqrt(Serr_val2[1] / S_val2[1]**2 + Serr_val2[0] / S_val2[0]**2) if S_val2[1] != 0 and S_val2[0] != 0 else 0


    result: Dict[str, Union[float, List[float]]] = {
        'Fptr': Fptr,
        'S_val': S_val,
        'Serr_val': Serr_val,
        'meanPdrift': meanPdrift,
        'meanTdrift': meanTdrift,
        'meanUdrift': meanUdrift,
        'S_val2': S_val2,
        'Serr_val2': Serr_val2,
        'Fha2o': Fha2o
    }

    return result



from typing import List, Dict, Union, Tuple

def ActrisProtocol_PICAB(path: str, name: str, ProtData: List[List[float]], masses: List[float], 
                         ProtEngData: List[List[float]], engnames: List[str], computerID: str) -> Dict[str, Union[float, List[float]]]:
    # Lines 18382-18386: Initialize variables
    helpi: int = 1
    T_FC: float = 31
    Fcal: float = 200
    
    # Lines 18388-18396: Adjust Fcal based on conditions
    if max(ProtEngData[0]) > 3770:
        Fcal = 240.0 * 298 / 273
    
    Fplus: float = 0
    if 'VOCUS PTR-TOF' in computerID:
        Fcal = 2000
    if len(computerID) == 0:
        Fcal = 2000
        Fplus = 5000 * 298 / 273  # Harvard
    
    # Lines 18398-18399: Calculate Flow
    Flow: float = float(Fcal) * 273 / 298 + Fplus + 9.6
    
    # Lines 18402-18403: Calculate Dil and T_loop
    Dil: float = (9.6 * 273 / (273 + T_FC)) / Flow
    T_loop: float = 273 + 100
    
    # Lines 18405-18411: Extract time and date information
    starttim: float = min(ProtEngData[:, 0])
    reltime: List[float] = [(t - starttim) * 24 * 60 for t in ProtEngData[:, 0]]
    PTRid: str = path.split('\\')[-2]
    Datumzeit: str = t09str(ProtEngData[0])  # Assuming t09str is defined elsewhere
    
    # Lines 18414-18432: Extract and process ion data
    I19 = ProtData[masses.index(19.017)]
    I21 = ProtData[masses.index(21.022)]
    I37 = ProtData[masses.index(37.028)]
    I38 = ProtData[masses.index(38.033)]
    I39 = ProtData[masses.index(39.033)]
    I34 = ProtData[masses.index(33.994)]
    I32 = ProtData[masses.index(31.989)]
    
    # Apply corrections to ion data
    I32 = [I34[i] * 242 if x > 6000 else x for i, x in enumerate(I32)]
    if max(I19) > 5000:
        I19 = [x * 487 for x in I21]
    I37 = [I38[i] * 645 if x < I38[i] * 645 else x for i, x in enumerate(I37)]
    
    # Lines 18436-18437: Calculate Prim
    Prim: List[float] = [i19 + i37 for i19, i37 in zip(I19, I37)]
    
    # Lines 18440-18441: Calculate F37 and F32
    F37: List[float] = [i37 / prim for i37, prim in zip(I37, Prim)]
    F32: List[float] = [i34 * 242 / prim for i34, prim in zip(I34, Prim)]
    
    # Lines 18443-18451: Extract drift data
    Pdrift: List[float] = ProtEngData[:, engnames.index('p-Drift')]
    H2Oscc: List[float] = ProtEngData[:, next(i for i, name in enumerate(engnames) if 'H2O scc' in name or 'H2O[scc' in name)]
    Udrift: List[float] = ProtEngData[:, engnames.index('Udrift')]
    Udx: List[float] = ProtEngData[:, engnames.index('Udx')]
    FCinlet: List[float] = ProtEngData[:, engnames.index('FCinlet')]
    Tdrift: List[float] = ProtEngData[:, engnames.index('T-Drift')]
    PC: List[float] = ProtEngData[:, next(i for i, name in enumerate(engnames) if 'PC mba' in name or 'PC[mba' in name)]
    
    # Lines 18453-18458: Adjust parameters for specific instruments
    if 'VOCUS PTR-TOF' in computerID:
        Pdrift = ProtEngData[:, engnames.index('SSQ pressure monitor [mbar]')]
        Tdrift = ProtEngData[:, engnames.index('TOF temp monitor')]
        Udrift = [x - 25 for x in ProtEngData[:, engnames.index('Spare 1 monitor [V]')]]
    elif len(computerID) == 0:  # PTR3
        Tdrift = ProtEngData[:, engnames.index('TPS body temp monitor')]
        Pdrift = ProtEngData[:, engnames.index('REAC   sensor: APR/CMR')]
    
    # Lines 18460-18524: Extract and process compound data
    D5_355 = ProtData[masses.index(355.07)] * 1.67
    D5_371 = ProtData[masses.index(371.102)] * 1.69
    D4_297 = ProtData[masses.index(297.083)] * 1.52
    D4_281 = ProtData[masses.index(281.052)] * 1.51
    D3_223 = ProtData[masses.index(223.064)] * 1.37
    D3_207 = ProtData[masses.index(207.033)] * 1.36
    benzCl3_181 = ProtData[masses.index(180.938)] * 2.48
    benzF3_133 = ProtData[masses.index(133.026)] * 1.07
    Pin_137 = ProtData[masses.index(137.132)] * 1.12
    Pin_81 = ProtData[masses.index(81.07)] * 1.07
    TMB_121 = ProtData[masses.index(121.101)] * 1.10
    XYL_107 = ProtData[masses.index(107.086)] * 1.09
    BEN_79 = ProtData[masses.index(79.054)] * 1.07
    MEK_73 = ProtData[masses.index(73.065)] * 1.05
    MVK_71 = ProtData[masses.index(71.049)] * 1.05
    MBO_87 = ProtData[masses.index(87.080)] * 1.06
    MBO_69 = ProtData[masses.index(69.070)] * 1.06
    MBO_41 = ProtData[masses.index(41.039)] * 1.03
    ACE_59 = ProtData[masses.index(59.049)] * 1.03
    ACETAL_45 = ProtData[masses.index(45.033)] * 1.02
    CH3CN_42 = ProtData[masses.index(42.034)] * 1.02
    METH_33 = ProtData[masses.index(33.033)] * 1.01
    BEN_78 = ProtData[masses.index(78.047)] * 1.07
    TMB_120 = ProtData[masses.index(120.094)] * 1.10
    XYL_106 = ProtData[masses.index(106.079)] * 1.09

    # Apply corrections for specific cases
    if abs(masses[masses.index(371.102)] - 371.102) > 1:
        D5_371 = D5_371 / 1000
    if abs(masses[masses.index(355.07)] - 355.07) > 1:
        D5_355 = D5_355 / 1000
    if abs(masses[masses.index(78.047)] - 78.047) > 1:
        BEN_78 = BEN_78 / 1000
    if abs(masses[masses.index(120.094)] - 120.094) > 1:
        TMB_120 = TMB_120 / 1000
    if abs(masses[masses.index(106.079)] - 106.079) > 1:
        XYL_106 = XYL_106 / 1000

    # Lines 18526-18533: Calculate TotSIG
    TotSIG: List[float] = [sum(x) for x in zip(METH_33, CH3CN_42, Acetal_45, Ace_59, MBO_87, MBO_69, MBO_41, MVK_71, MEK_73, Ben_79,
                                               Xyl_107, TMB_121, Pin_137, Pin_81, benzF3_133, benzCl3_181, D3_223, D3_207, D4_297, D4_281,
                                               D5_371, D5_355)]
        
    # Lines 18536-18559: Extract isotope data
    D5_358iso = ProtData[masses.index(358.073)] / 0.0358
    D5_374iso = ProtData[masses.index(374.105)] / 0.0370
    D4_299iso = ProtData[masses.index(299.086)] / 0.1178
    D4_283iso = ProtData[masses.index(283.055)] / 0.1170
    D3_225iso = ProtData[masses.index(225.067)] / 0.0923
    D3_209iso = ProtData[masses.index(209.036)] / 0.0916
    benzCl3_182iso = ProtData[masses.index(181.941)] / 0.0267
    benzF3_134iso = ProtData[masses.index(134.029)] / 0.0612
    Pin_138iso = ProtData[masses.index(138.135)] / 0.0987
    Pin_82iso = ProtData[masses.index(82.073)] / 0.0617
    TMB_122iso = ProtData[masses.index(122.104)] / 0.0896
    XYL_108iso = ProtData[masses.index(108.089)] / 0.0804
    BEN_80iso = ProtData[masses.index(80.057)] / 0.0615
    MEK_74iso = ProtData[masses.index(74.068)] / 0.0426
    MVK_72iso = ProtData[masses.index(72.052)] / 0.0424
    MBO_88iso = ProtData[masses.index(88.083)] / 0.0526
    MBO_70iso = ProtData[masses.index(70.073)] / 0.0522
    MBO_42iso = ProtData[masses.index(42.042)] / 0.0320
    ACE_60iso = ProtData[masses.index(60.052)] / 0.0325
    ACETAL_46iso = ProtData[masses.index(46.036)] / 0.0220
    CH3CN_43iso = ProtData[masses.index(43.037)] / 0.0251
    METH_34iso = ProtData[masses.index(34.036)] / 0.0116

    MBO_42iso *= 0.000001
    METH_34iso *= 0.00001
     
    # Lines 18562-18567: Calculate TotSIGiso
    TotSIGiso: List[float] = [sum(x) for x in zip(METH_34iso, CH3CN_43iso, Acetal_46iso, Ace_60iso, MBO_88iso, MBO_70iso, MBO_42iso, MVK_72iso, MEK_74iso, Ben_80iso,
                                                  Xyl_108iso, TMB_122iso, Pin_138iso, Pin_82iso, benzF3_134iso, benzCl3_182iso, D3_225iso, D3_209iso, D4_299iso, D4_283iso,
                                                  D5_374iso, D5_358iso)]
    
    # Lines 18570-18574: Define mTOz and compounds lists
    mTOz: List[float] = [33, 42, 45, 59, 87, 69, 41, 71, 73, 79, 107, 121, 137, 81, 133, 181, 223, 207, 297, 281, 371, 355, -1, -1, -1, -1, -1, -1, 19, 37, 39, 34, 43, 46, 60, 88, 70, 42, 72, 74, 80, 108, 122, 138, 82, 134, 182, 225, 209, 299, 283, 374, 358, -1, -1, -1, 78, 32, 120, 106]
    compounds: List[str] = ['METH_33', 'CH3CN_42', 'Acetal_45', 'Ace_59', 'MBO_87', 'ISO/MBO_69', 'ISO/MBO_41', 'MVK_71', 'MEK_73', 'Ben_79',
                            'Xyl_107', 'TMB_121', 'MT_137', 'MT_81', 'benzF3_133', 'benzCl3_181', 'D3_223', 'D3_207', 'D4_297', 'D4_281',
                            'D5_371', 'D5_355', 'TotSIG', 'FCinlet', 'reltime', 'TotSIGiso', 'prim', 'Tdrift', 'I19', 'I37',
                            'I39', 'METH_34iso', 'CH3CN_43iso', 'Acetal_46iso', 'Ace_60iso', 'MBO_88iso', 'MBO_70iso', 'MBO_42iso', 'MVK_72iso', 'MEK_74iso', 'Ben_80iso',
                            'Xyl_108iso', 'TMB_122iso', 'Pin_138iso', 'Pin_82iso', 'benzF3_134iso', 'benzCl3_182iso', 'D3_225iso', 'D3_209iso', 'D4_299iso', 'D4_283iso',
                            'D5_374iso', 'D5_358iso', 'Udx', 'pdrift', 'Udrift', 'BEN_78', 'I32', 'TMB_120', 'XYL_106']
    
    # Lines 18576-18579: Define index variables
    i_Tdrift: int = 27
    i_Udx: int = 53
    i_pdrift: int = 54
    i_Udrift: int = 55
    
    # Lines 18586-18920: Injection detection and processing
    LongDUR = 10
    tims = [0]
    timmax = [0]
    timend = [0]
    sigmax = [0]
    SearchStart = 1
    SearchEnd = 0

    data = ACE_59
    dataSTR = 'ACE_59'
    if computerID == '2197a787-6c47-4c78-88d1-28fb6598bcc2' and abs(min(ProtEngData[:, 0]) - 3184.484) < 0.2:
        data = MEK_73
        dataSTR = 'MEK_73'
    if computerID == '50b4c522-8345-4094-bb5d-683484b5253e' and abs(min(ProtEngData[:, 0]) - 3190.6) < 0.2:
        data = benzF3_133
        dataSTR = 'benzF3_133'
    if len(computerID) == 0:
        data = CH3CN_42
        dataSTR = 'CH3CN_42'

    mini = np.quantile(data, 0.1)
    midi = np.quantile(data, 0.4)
    i30s = int(0.5 * len(data) / max(reltime))
    faktor2 = 3 if computerID == '2197a787-6c47-4c78-88d1-28fb6598bcc2' else 2
    step = 5 if len(computerID) == 0 else 1

    for i in range(3, len(data) - 2):
        if SearchStart == 1:
            if abs(data[i] - mini) / abs(data[i - 2 * step] - mini) > 3 and data[i] > faktor2 * midi:
                tims.append(i)
                istart = i
                mini = min(data[i - step], data[i - 2 * step], data[i - 3 * step])
                i += 1 + step - 1
                SearchStart = 0
                SearchEnd = 1

        if i > i30s:
            mininew = min(data[i - i30s:i])
            if mininew < midi * faktor:
                mini = mininew

        if SearchEnd == 1:
            if abs(data[i - step] - mini) / abs(data[i] - mini) > 1.5 and data[i - step] > midi:
                timend.append(i)
                sigmax.append(max(data[istart:i]) / midi)
                helpi = data[istart:i].index(max(data[istart:i]))
                timmax.append(istart + helpi)
                i += 1
                istart += 1
                SearchStart = 1
                SearchEnd = 0

    tims = tims[1:]
    timend = timend[1:]
    sigmax = sigmax[1:]
    timmax = timmax[1:]

    # Clean and process detected injections
    injectiontimes = [reltime[t] for t in tims]
    injectiontimesEND = [reltime[t] for t in timend]
    duration = [(end - start) * 60 for start, end in zip(injectiontimes, injectiontimesEND)]

    # Clean short injections
    Shortsig = np.quantile([s for s, d in zip(sigmax, duration) if d < LongDUR], 0.6)
    intex = [i for i, (s, d) in enumerate(zip(sigmax, duration)) if d < LongDUR and s < Shortsig / 5]
    for i in sorted(intex, reverse=True):
        del tims[i], timend[i], injectiontimes[i], injectiontimesEND[i], duration[i], sigmax[i], timmax[i]

    # Process sequences
    PREVdistance = [-9999] + [-60.0 * (injectiontimes[i] - injectiontimes[i+1]) for i in range(len(injectiontimes)-1)]
    NEXTdistance = [60.0 * (injectiontimes[i+1] - injectiontimes[i]) for i in range(len(injectiontimes)-1)] + [-9999]
    sequence = [100 if duration[0] > LongDUR else 1]
    for i in range(1, len(tims)):
        if duration[i] <= LongDUR and abs(PREVdistance[i]) < 8.5 and sequence[i-1] < 99:
            sequence.append(sequence[i-1] + 1)
        elif duration[i] > LongDUR:
            sequence.append(100)
        else:
            sequence.append(1)

    for i in range(len(tims)):
        if sequence[i] == 10:
            for j in range(10):
                sequence[i-j] = 10

    # Clean long injections
    Longsig = max([s for s, seq in zip(sigmax, sequence) if seq == 100])
    intex = [i for i, (s, seq) in enumerate(zip(sigmax, sequence)) if seq == 100 and s < Longsig / 10]
    for i in intex:
        sequence[i] = -9999

    # Final cleaning
    if max([seq for seq in sequence if seq not in [10, 100]]) > -0.5:
        valid_indices = [i for i, seq in enumerate(sequence) if seq > -9999]
        tims = [tims[i] for i in valid_indices]
        timend = [timend[i] for i in valid_indices]
        injectiontimes = [injectiontimes[i] for i in valid_indices]
        injectiontimesEND = [injectiontimesEND[i] for i in valid_indices]
        duration = [duration[i] for i in valid_indices]
        sigmax = [sigmax[i] for i in valid_indices]
        timmax = [timmax[i] for i in valid_indices]
        sequence = [seq for seq in sequence if seq > -9999]

    # Separate short and long injections
    injectSHORT = [t for t, d in zip(injectiontimes, duration) if d < 4]
    injectSHORTend = [t for t, d in zip(injectiontimesEND, duration) if d < 4]
    timsSHORT = [t for t, d in zip(tims, duration) if d < 4]
    timsSHORTend = [t for t, d in zip(timend, duration) if d < 4]

    injectLONG = [t for t, d in zip(injectiontimes, duration) if d > 4]
    injectLONGend = [t for t, d in zip(injectiontimesEND, duration) if d > 4]
    timsLONG = [t for t, d in zip(tims, duration) if d > 4]
    timsLONGend = [t for t, d in zip(timend, duration) if d > 4]


    # Lines 18922-19920: Additional processing and calculations
    # Calculate sensitivities
    sensitivities = calculate_sensitivities(RRK, t_react, p_d, T_d)

    # Transmission retrieval
    bestEX, MASSlow, WIDTHlow, MASShigh, WIDTHhigh = retrieve_transmission(counts, countsBG)

    # Calculate normalized counts
    NC = calcNC(bestEX, MASSlow, WIDTHlow, MASShigh, WIDTHhigh, counts, countsBG)

    # Calculate essentials
    essen = calculate_essentials(NC, Fmol, en)

    # Plot transmission curve
    plot_transmission(bestEX, MASSlow, WIDTHlow, MASShigh, WIDTHhigh, essen, sensitivities)

    # Calculate and plot XR data
    XRdata = calculate_XR(countsBGTR, countsTR)
    plot_XR_data(XRdata)

    # Process short injections
    process_short_injections(CompCUB, timestep)

    # Export data
    export_sens_data(path2, PTRid, Datumzeit, KINSPECS, ncpsPERppbTR, ncpsPERppbTRxr, cpsPERppbRAW)
    export_XR_data(path2, PTRid, Datumzeit, masslow, bestEX, widthlow, masshigh, widthhigh, XRdata)

    # Generate final plots
    plot_basic_parameters(reltime, TotSIG, I19, I37, data)
    plot_fast_injections(CompCUB, VMRppm, dil)

    # Create and return result structure (simplified for this example)
    result: Dict[str, Union[float, List[float]]] = {
        'Fptr': [1, 1, 1],
        'S_nam': ['time'] + compounds[:22],
        'S_val': [starttim] + S_val[:22],
        'Serr_val': [0] + Serr_val[:22]
    }
    return result

def ActrisProtocol_nov_2018():
    #todo
    pass

def ActrisProtocol_oct_2018():
    ### TODO
    pass

def ActrisProtocol_03_2018():
    ### TODO
    pass
