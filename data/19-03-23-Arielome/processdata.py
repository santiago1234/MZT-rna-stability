import os
import pandas as pd
from mztrnastabilitty.dataprocessing import arielome


def remove_redundata_cols(processed_file):
    return processed_file.drop(['sequence', 'stop_position'], axis=1)


def processing(rawarielome):
    '''
    process arielome data
    '''
    par = arielome.processarielome(rawarielome, remove_redundata_cols)
    # I will split sequences with premature and not premature stop codon
    with_premature_stop = par[par.utr3.str.len() != 0]
    no_stop_codon = par[par.utr3.str.len() == 0]
    no_stop_codon = no_stop_codon.drop('utr3', axis=1)
    return with_premature_stop, no_stop_codon


def write_file(data, outpath):
    """Write file, append if exists"""
    if not os.path.isfile(outpath):
        data.to_csv(outpath, index=False)
    else:
        data.to_csv(outpath, mode='a', header=False, index=False)


def processfile(fileid):
    """
    process data
    """
    print('processing file: {}'.format(fileid))
    filechunks = arielome.getrawfile(filepaths[fileid])
    # create out dir
    outdirname = 'processed_arielome_data'
    if not os.path.exists(outdirname):
        os.makedirs(outdirname)
    outname = [os.path.join(outdirname, (fileid + '_' + x + '.csv'))
               for x in ('premature_stop', 'no_stop')]
    for dchunk in filechunks:
        with_premature_stop, no_stop_codon = processing(dchunk)
        # write files
        write_file(with_premature_stop, outname[0])
        write_file(no_stop_codon, outname[1])
    return None


# run the processing
filepaths = arielome.arielome_dpaths("../../")
for fileid in filepaths:
    processfile(fileid)
