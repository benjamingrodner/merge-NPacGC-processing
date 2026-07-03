import os
import gzip


def main():
    dir_expcounts = '/mnt/nfs/projects/armbrust-metat/gradients2/g2_dcm_rr_pa_metat/kallisto'
    dict_batch_fn = {
        'G2PA.RR':f'{dir_expcounts}/G2Inc_mappedReadAbundance.tsv.gz',
        'G2PA.DCM':f'{dir_expcounts}/G2DCM_mappedReadAbundance.tsv.gz',
    }
    dir_splitexp = '/scratch/bgrodner/metat_data'
    header = 'target_id\tlength\teff_length\test_counts\ttpm|name_file\n'
    dict_sn_handle = {}
    for batch, fn in dict_batch_fn.items():
        dir_out = f'{dir_splitexp}/kallisto_counts/{batch}'
        os.makedirs(dir_out, exist_ok=True)
        dict_sn_fn = {}
        print('Splitting file:', batch)
        i = 0
        with gzip.open(fn, 'rt') as f:
            for line in f:
                l = line.strip().split('\t')
                if (l[0] != 'target_id') and (len(l[0]) > 0):
                    try:
                        sn = l[4].split('|')[1]
                    except:
                        raise IndexError(f'Line = {line}\nl = {l}')
                    fn_out = f'{dir_out}/{batch}-kallisto_counts-{sn}'
                    handle = dict_sn_handle.get(sn)
                    if handle is None:
                        handle = open(fn_out, 'w')
                        handle.write(header)
                        dict_sn_handle[sn] = handle
                    else:
                        handle = dict_sn_handle[sn]
                    handle.write(line)
                i += 1
                if i % 1e6 == 0:
                    print(f'Lines read: {i}', end='\r')

    for sn, handle in dict_sn_handle.items():
        handle.close()

    return

if __name__ == "__main__":
    main()
