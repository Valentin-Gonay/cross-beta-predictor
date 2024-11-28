''' 
Valentin Gonay
17/10/23
'''


# Analyze sequence and create / extract features from it

# Extract features from sequence #
def get_features(
        sequence: str, 
        feature_list: list, 
        classification: str
        ):
    '''Create the wanted feature of the given protein sequence for prediction
    
    :param sequence: The protein sequence you want to predict
    :type sequence: str 

    :param feature_list: The list of feature to create for the relative sequence
    :type feature_list: list 

    :param classification: The classification method for converting AA into a group value. 
    Can be: 'mode_1', 'mode_2' or 'mode_3'
    :type classification: str 

    :return: The created entry, ready to be predict
    :rtype: list
    '''

    grp_seq = convert_AA_grp(sequence, classification)

    # get AA compo
    AA_list = get_AA_list(feature_list)
    AA_compo = get_AA_compo(sequence, AA_list)

    # get groups transi
    transi_list = get_transi_list(feature_list)
    transi_compo = get_transi_compo(grp_seq, transi_list) 

    # get groups compo
    grp_list = get_grp_list(feature_list)
    grp_compo = get_grp_compo(grp_seq, grp_list)
    
    return AA_compo + transi_compo + grp_compo



# Analyse prediction result by amino acids #
# get AR
def find_amyloid_region(
        confidence_mean_list: list, 
        threshold: float = 0.5
        ):
    '''Get all predicted amyloid regions from a confidence list using the given threshold
    
    :param confidence_mean_list: The list of confidence mean for each amino acids
    :type confidence_mean_list: list 

    :param threshold: The threshold to use for the confidence score, values higher than the 
    treshold will be considered positive. If not provided, defaults to 0.5
    :type threshold: float 

    :return: The list of amyloidogenic regions (= None if no regions are detected). Note that 
    regions has to have a length superior or equal to 15 amino acids
    :rtype: list
    '''

    is_AR_start = False
    AR_start_pos = None
    region_list = []
    for i in range(len(confidence_mean_list)):
        current_score = list(confidence_mean_list[i].values())[0]
        if i == len(confidence_mean_list)-1:
            next_score = None
        else:
            next_score = list(confidence_mean_list[i+1].values())[0]
        
        if not is_AR_start:
            is_AR_start = find_AR_start(current_score, threshold)
            if is_AR_start:
                AR_start_pos = i
        if is_AR_start:
            is_AR_end = find_AR_end(current_score, next_score, threshold)
            if is_AR_end:
                AR_end_pos = i
                if AR_end_pos - AR_start_pos >= 14:
                    # +1 in order to fit with real sequence numerotation (start at 1 and not 0)
                    region_list.append([AR_start_pos+1,AR_end_pos+1]) 
                    is_AR_start = False
                    is_AR_end = False
                else:
                    is_AR_start = False
                    is_AR_end = False
    return region_list



def find_AR_start(
        AA_score: float, 
        threshold: float
        ):
    '''Check the score of a amino acid to see if it can be the beggining of a new AR
    
    :param AA_score: The score of the current amino acid
    :type AA_score: float 

    :param threshold: The threshold to use for the confidence score, 
    values higher than the treshold will be considered positive
    :type threshold: float 
    
    :return: True if it can be the beggining of a AR, else False
    :rtype: bool
    '''

    return AA_score > threshold



def find_AR_end(
        AA_score: float, 
        next_AA_score: float, 
        threshold: float
        ):
    '''Check the score of a amino acid to see if it can be the end of a new AR
    
    :param AA_score: The score of the current amino acid.
    :type AA_score: float 

    :param next_AA_score: The score of the next amino acid.
    :type next_AA_score: float 

    :param threshold: The threshold to use for the confidence score, values higher than the treshold 
    will be considered positive
    :type threshold: float
    
    :return: True if it can be the end of a AR, else False
    :rtype: bool
    '''

    if next_AA_score == None:
        return AA_score > threshold
    return AA_score > threshold and next_AA_score <= threshold



# Find optimum window length for the prediction by amino acids #
def get_length_threshold(sequence: str):
    '''Get the default region length for prediction of AR in long sequence ( > 150 amino acids)
    
    :param sequence: The sequence use for the prediction
    :type sequence: str 
    
    :return: The length of the AR for the prediction 
    :rtype: int
    '''

    length_seq = len(sequence)

    if length_seq <=150:
        return 15
    elif length_seq >= 500:
        return 50
    else:
        return int(length_seq/10)



# Compute the Amino acids composition of one sequence #
def get_AA_list(feature_list: list):
    '''Extract the list of Amino Acid present in the given feature list
    
    :param feature_list: The given feature list
    :type feature_list: list 
    
    :return: The list of Amino Acids present in the feature list
    :rtype: list
    '''

    res = []
    for i in range(len(feature_list)):
        if len(feature_list[i]) == 1:
            res.append(feature_list[i])
    return res



def get_AA_compo(
        sequence: str, 
        AA_list: list
        ):
    '''Get the Amino Acid (AA) composition of the sequence, only compute the composition of the AA 
    present in the given list
    
    :param sequence: The protein sequence
    :type sequence: str 

    :param AA_list: The list of amino acids to get the composition
    :type AA_list: list 
    
    :return: The composition of amino acid (in the same order than the given AA_list)
    :rtype: list
    '''

    AA_compo = []
    for AA in AA_list:
        AA_compo.append(get_nb_of_1_elem(sequence,AA))
    return AA_compo



def get_nb_of_1_elem(
        seq: str,
        elem: str
        ):
    '''Count the number of each given element and give their proportion
    
    :param seq: The AA sequence
    :type seq: str 

    :param elem: The element to count
    :type elem: str 

    :return: The proportion of the given AA in the given sequence
    :rtype: float
    '''

    return round(seq.count(elem)/len(seq),4)



# Create Group sequence based on a classification method #
def convert_AA_grp(
        sequence: str, 
        classification: str
        ):
    '''Convert the Amino acids (AA) of the given sequence into their groups according to the 
    given classification method
    
    :param sequence: The protein sequence
    :type sequence: str 

    :param classification: The classification method for converting AA into a group value. 
    Can be: 'mode_1', 'mode_2' or 'mode_3'
    :type classification: str 
    
    :return: The sequence of groups value
    :rtype: str
    '''

    if classification == 'mode_1':
        return get_CandF_seq(sequence)
    elif classification == 'mode_2':
        return get_group_seq(sequence)
    elif classification == 'mode_3':
        return get_group_seq_2(sequence)
    else:
        print("\nError, the given classification method is not recognized by the script.",
              "Please retry with one of the following method:\n'mode_1', 'mode_2', or 'mode_3'\n")
        return None
    


# Classification methods #
# C_and_F
def get_CandF_seq(seq: str):
    '''Create a sequence by replacing the amino acids of the given sequence by its group belonging
    
    :param seq: The amino acid sequence
    :type seq: str 
    
    :return: The AA group sequence
    :rtype: str
    '''

    A_group = ['R', 'K', 'E', 'D', 'Q', 'N', 'B', 'Z'] # polar
    B_group = ['G', 'A', 'S', 'T', 'P', 'H', 'Y'] # neutral
    C_group = ['C', 'V', 'L', 'I', 'M', 'F', 'W', 'J'] # hydrophobic
    X_group = ['X', 'O', 'U'] # other

    converted_seq = ''
    for aa in seq:
        if aa in A_group:
            converted_seq += 'A'
        elif aa in B_group:
            converted_seq += 'B'
        elif aa in C_group:
            converted_seq += 'C'
        elif aa in X_group:
            converted_seq += 'X'
    return converted_seq 



# Perso
def get_group_seq(seq: str):
    '''Create a sequence by replacing the amino acids of the given sequence by its group belonging
    
    :param seq: The amino acid sequence
    :type seq: str 
    
    :return: The A.A. group sequence
    :rtype: str
    '''

    P_group = ['P']
    G_group = ['G']
    A_group = ['H', 'N', 'T', 'Q', 'N', 'C', 'S'] # polar
    B_group = ['K', 'R', 'E', 'D'] # charged
    C_group = ['I', 'L', 'M', 'V', 'W', 'Y', 'F', 'A'] # hydrophobic
    X_group = ['X', 'O', 'U', 'B', 'Z', 'J'] # other

    converted_seq = ''
    for aa in seq:
        if aa in A_group:
            converted_seq += 'A'
        elif aa in B_group:
            converted_seq += 'B'
        elif aa in C_group:
            converted_seq += 'C'
        elif aa in X_group:
            converted_seq += 'X'
        elif aa in P_group:
            converted_seq += 'P'
        elif aa in G_group:
            converted_seq += 'G'
    return converted_seq



# Perso_2
def get_group_seq_2(seq: str):
    '''Create a sequence by replacing the amino acids of the given sequence by its group belonging
    
    :param seq: The amino acid sequence
    :type seq: str 
    
    :return: The A.A. group sequence
    :rtype: str
    '''

    P_group = ['P']
    G_group = ['G']
    A_group = ['H', 'T', 'C', 'S'] # polar
    B_group = ['K', 'R', 'E', 'D'] # charged
    C_group = ['I', 'L', 'M', 'V', 'W', 'Y', 'F', 'A'] # hydrophobic
    D_group = ['Q', 'N'] # Important for amyloid formation
    X_group = ['X', 'O', 'U', 'B', 'Z', 'J'] # other

    converted_seq = ''
    for aa in seq:
        if aa in A_group:
            converted_seq += 'A'
        elif aa in B_group:
            converted_seq += 'B'
        elif aa in C_group:
            converted_seq += 'C'
        elif aa in X_group:
            converted_seq += 'X'
        elif aa in P_group:
            converted_seq += 'P'
        elif aa in G_group:
            converted_seq += 'G'
        elif aa in D_group:
            converted_seq += 'D'
    
    return converted_seq



# Compute the group composition of one sequence #
def get_grp_list(feature_list: list):
    '''Extract the list of group present in the given feature list
    
    :param feature_list: The given feature list
    :type feature_list: list 
    
    :return: The list of group present in the feature list
    :rtype: list
    '''

    res = []
    for i in range(len(feature_list)):
        if 'grp_' in feature_list[i]:
            res.append(feature_list[i])
    return res



def get_grp_compo(
        sequence: str, 
        grp_list: list
        ):
    '''Get the group composition of the sequence, only compute the composition of the group 
    present in the given list
    
    :param sequence: The protein sequence (with amino acid groups)
    :type sequence: str 

    :param grp_list: The list of group to get the composition
    :type grp_list: list 
    
    :return: The composition of Amino acid groups (in the same order than the given grp_list)
    :rtype: list
    '''

    grp_compo = []
    for grp in grp_list:
        grp_compo.append(get_nb_of_1_elem(sequence,grp))
    return grp_compo



# Compute the transition between each groups in one sequence #
def get_transi_list(feature_list: list):
    '''Extract the list of group transition present in the given feature list
    
    :param feature_list: The given feature list
    :type feature_list: list 

    :return: The list of group transition present in the feature list
    :rtype: list
    '''
    
    res = []
    for i in range(len(feature_list)):
        if '_to_' in feature_list[i]:
            res.append(feature_list[i])
    return res


def get_transi_compo(
        sequence: str, 
        transi_list: list
        ):
    '''Get the group composition of the sequence, only compute the composition of the transition 
    between group present in the given list
    
    :param sequence: Description
    :type sequence: str 

    :param transi_list: Description
    :type transi_list: list

    :return: The composition of transition between groups (in the same order than 
    the given transi_list)
    :rtype: list
    '''

    transi_compo = []
    for transi in transi_list:
        transi_compo.append(get_one_transi_compo(sequence,transi))
    return transi_compo


def get_one_transi_compo(
        sequence: str, 
        transi: str
        ):
    '''Compute the composition of one given transition for the given sequence
    
    :param sequence: The protein sequence (with amino acid groups)
    :type sequence: str 

    :param transi: The specific transition to count and extract the composition 
    (Must be in format 'X_to_X')
    :type transi: str 
    
    :return: The composition of the given transition in the sequence
    :rtype: list
    '''
    
    parse_transi = transi.split('_to_')
    sequence_length = len(sequence)
    count = 0
    AA_1 = parse_transi[0]
    AA_2 = parse_transi[1]
    for i in range(len(sequence)-1):
        if sequence[i] == AA_1 and sequence[i+1] == AA_2:
            count += 1
    return count/sequence_length

