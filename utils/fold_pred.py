''' 
Valentin Gonay
16/10/23
'''

# INTERN IMPORTS
import utils.iupred3.iupred3_lib as iupred

# EXTERN IMPORTS
import numpy as np







### Compute structurality prediction (IUPred3) ###
# IUPred prediction (disorder score)
def get_iupred_res_df(
        dataframe, 
        sequence_col: str, 
        ar_col : str | None = None
        ):
    '''Get the iupred score of the given sequence and store it in the proper column of 
    the given dataframe.
    
    :param dataframe: The dataframe where to extract the sequence and store the result
    :type dataframe: DataFrame 

    :param sequence_col: The name of the column containing the sequence in the given dataframe
    :type sequence_col: str 

    :param ar_col: The column where to find the amyloid region of the given protein. If given, the \
        IUPred score will be the mean score of this region. If None, the score will be the mean of the \
        all protein (Default = None)
    :type ar_col: str | None

    :return: The updated dataframe with the new IUPred value
    :rtype: DataFrame
    '''

    # Initialize start time
    for i, row in dataframe.iterrows():
        seq = row[sequence_col]
        iupred_res = iupred.iupred(seq, "long")
        if ar_col != None:
            ar = row[ar_col]
            ar_split = ar.split('-')
            mean_res = np.mean(iupred_res[0][int(ar_split[0]):int(ar_split[1])+1])
        else:
            mean_res = np.mean(iupred_res[0])
        dataframe.at[i, 'IUPred3_prediction'] = mean_res
    return dataframe


# Global fold score function
def get_IUPred_allprot(sequence: str):
    '''Get the list of value comming from IUPred for the given protein sequence.
    
    :param sequence: The protein sequence that is only composed of the 20 essential amino acids code
    :type sequence: str 
    
    :return: The dictionary with the lists of value produced by IUPred (key = 'iupred_score_list')
    :rtype: dict
    '''
    
    if len(sequence) < 19:
        iupred_res = iupred.iupred(sequence, "short", smoothing='strong')[0]
    else:
        iupred_res = iupred.iupred(sequence, "long")[0]
    return {
        'iupred_score_list':iupred_res,
    }
