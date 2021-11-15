import itertools
import pandas as pd
import numpy as np
import seaborn as sns
import matplotlib.pyplot as plt
import scanpy as sc
import episcanpy as epi

def top_markers(adata, ntop=100):
    
    markers = None

    # If processed with Scanpy
    if 'rank_genes_groups' in adata.uns:
        markers = pd.DataFrame(adata.uns['rank_genes_groups']['names']).head(ntop)
        
    # Else if processed with Episcanpy
    elif 'rank_features_groups' in adata.uns:
        markers = pd.DataFrame(adata.uns['rank_features_groups']['names']).head(ntop)
    else:
        print('No rank_genes_groups variable found in the submitted anndata.')
        print('Please run scanpy.tl.rank_genes_groups to find group markers.')
        
    return(markers)



def plot_jaccardindex(data, xlab='x', ylab='y'):
    
    data = data.round(decimals=2)
    # This will use the same order for the labels
    data = data.sort_index(0).sort_index(1)
    # This will reverse the matrix so the diagonal looks good
    data = data[::-1]
    
    plot = sns.heatmap(data, annot = True, vmin = 0., vmax = 1., cmap = 'Reds')
    plt.xticks(rotation = 70)
    plt.xlabel(xlab, fontsize = 30)
    plt.ylabel(ylab, fontsize = 30)
    plt.title('Jaccard Index scores', fontsize = 40, pad = 20)
    plt.show()
    
    return(plot)


def matchscore(ref_markers, obs_markers, plot=True):
    
    score = 0
    mat = []
    
    for ref_cluster in ref_markers:
        ref_len = len(ref_markers[ref_cluster])
        
        jacc_index = []
        
        for obs_cluster in obs_markers:
            
            obs_len = len(obs_markers[obs_cluster])
            shared = len(set(ref_markers[ref_cluster]).intersection(obs_markers[obs_cluster]))
            jacc = shared / (ref_len + obs_len - shared)
            jacc_index.append(jacc)
            
        score += np.max(jacc_index)
        mat.append(jacc_index)
    
    mat = pd.DataFrame(mat, index = list(ref_markers.columns),
                       columns = list(obs_markers.columns))
    
    score = score / ref_markers.shape[1]
    max_ji = mat.max()
    anno_lab = mat.idxmax()
    if plot:
        plot_jaccardindex(mat, xlab='My dataset', ylab='Reference dataset')
    
    return((score, anno_lab, max_ji, mat))




def find_common_genes(ref, sample, target_n_genes=2000):
    
    print(f'Dimensions before filtering: {ref.shape} and {sample.shape}\n')
    
    n = 5000
    common = []
    while len(common) < target_n_genes:
        
        print(f'Looking for {n} HVG')
        sc.pp.highly_variable_genes(ref, n_top_genes=n)
        sc.pp.highly_variable_genes(sample, n_top_genes=n)

        ref_genes = ref.var.highly_variable[ref.var.highly_variable == True].index
        sample_genes = sample.var.highly_variable[sample.var.highly_variable == True].index

        common = list(set(sample_genes).intersection(ref_genes))
        print(f'Found {len(common)} genes in common')
        n += 1000

    print('Filtering the data to these features and scaling')
    
    ref = ref[:, common]
    sample = sample[:, common]
    sc.pp.scale(ref, max_value=10)
    sc.pp.scale(sample, max_value=10)
    print(f'\nDimensions after filtering: {ref.shape} and {sample.shape}')
    
    return ref, sample



def find_common_variable_peaks(ref, sample, target_n_peaks=5000):
    
    print(f'Dimensions before filtering: {ref.shape} and {sample.shape}\n')
    
    n = target_n_peaks
    common = []
    while len(common) < target_n_peaks:
        
        print(f'Looking for {n} HVP')
        ref_var = epi.pp.select_var_feature(ref, nb_features=n, show=False, copy=True)
        sample_var = epi.pp.select_var_feature(sample, nb_features=n, show=False, copy=True)

        ref_genes = list(ref_var.var.index)
        sample_genes = list(sample_var.var.index)

        common = list(set(sample_genes).intersection(ref_genes))
        print(f'Found {len(common)} peaks in common')
        n += 1000

    print('Filtering the data to these features and scaling')
    
    ref = ref[:, common]
    sample = sample[:, common]
    sc.pp.scale(ref, max_value=10)
    sc.pp.scale(sample, max_value=10)
    print(f'\nDimensions after filtering: {ref.shape} and {sample.shape}')
    
    return ref, sample



def find_refmarkers_in_variable_peaks(ref, sample, target_n_peaks=5000):
    
    print(f'Dimensions before filtering: {ref.shape} and {sample.shape}\n')
    
    dfn = ref.obs['predicted.id'].value_counts() > 20
    keep_celltypes = dfn[dfn == True].index.to_list()
    ref = ref[ref.obs['predicted.id'].isin(keep_celltypes)]
    sc.pp.normalize_total(ref)
    epi.pp.log1p(ref)
    
    sc.pp.normalize_total(sample)
    epi.pp.log1p(sample)
    
    n_cat = len(ref.obs['predicted.id'].cat.categories)
    n = int(round(target_n_peaks*2/n_cat))
    common = []
    
    while len(common) < target_n_peaks:
        
        print(f'Looking for {n} marker features from each label')
        epi.tl.rank_features(ref, 'predicted.id', omic="ATAC", use_raw=False, n_features=n)
        marker_peaks = list(ref.uns['rank_features_groups']['names'])
        marker_peaks = list(itertools.chain.from_iterable(marker_peaks))
        print(len(marker_peaks))

        sample_peaks = epi.pp.select_var_feature(sample, nb_features=len(marker_peaks), show=False, copy=True)
        sample_peaks = list(sample_peaks.var.index)

        common = list(set(marker_peaks).intersection(sample_peaks))
        print(f'Found {len(common)} peaks in common')
        n += 100

    print('Filtering the data to these features and scaling')
    
    ref = ref[:, common]
    sample = sample[:, common]
    sc.pp.scale(ref, max_value=10)
    sc.pp.scale(sample, max_value=10)
    print(f'\nDimensions after filtering: {ref.shape} and {sample.shape}')
    
    return ref, sample



def filter_peaks_width(peakdic, width_range):

    minim, maxim = width_range
    before, after = 0, 0
    peak_lengths = []
    
    for seq in peakdic:
        before += len(peakdic[seq])
        
        for peak in peakdic[seq]:
            peak_lengths.append(len(peak))
            
            if len(peak) < minim or len(peak) > maxim:
                peakdic[seq].remove(peak)
                
        after += len(peakdic[seq])
    
    plt.hist(np.log10(peak_lengths), bins = 100)
    plt.axvline(np.log10(minim), color='k', 
                linestyle='dashed', linewidth=1)
    plt.axvline(np.log10(maxim), color='k', 
                linestyle='dashed', linewidth=1)
    plt.show()
    
    print(f'{before} peaks before filtering')
    print(f'{after} peaks after filtering')
    return(peakdic)



def peak2range(peak):
    
    chrom, r = peak.split(':')
    r = r.split('-')
    rang = range(int(r[0]), int(r[1])+1)
    return chrom, rang
    
    
    
def range2peak(chrom, rang):
    return f'{chrom}:{rang.start}-{rang.stop-1}'
    
    
    
def generate_peakdic(peak_list):
    
    dic = {}
    for peak in peak_list:

        chrom, rang = peak2range(peak)
        
        if chrom in dic:
            dic[chrom].append(rang)
        else:
            dic[chrom] = [rang]

    return dic
    
    
    
def merge_ranges(peak_list):
    
    start, stop = None, None
    
    for p in peak_list:
        # First peak in the list
        if start is None:
            start, stop = p.start, p.stop

        # If previous peak end is less than new peak start, 
        # they're distinct peaks. Yield and save new:
        elif stop < p.start:
            yield range(start, stop)
            start, stop = p.start, p.stop

        # If the stop from previous peak is after new peak start,
        # they should be the same peak. Save for posterior merge
        else:
            stop = p.stop
        
    yield range(start, stop)

    

def filter_chroms(chroms):
    return [x for x in chroms if x.startswith('chr')]
    

    
def get_common_peaks(adatas, only_standard_chroms=False):
    
    peaks_dic_list = [generate_peakdic(a.var.index.tolist())
                      for a in adatas]
    
    common_seqs = set.intersection(*map(set, peaks_dic_list))
    if only_standard_chroms:
        common_seqs = filter_chroms(common_seqs)
 
    common_peaks = {}
    
    for s in common_seqs:
        
        seq_lens = [len(adata[s]) for adata in peaks_dic_list]
        if min(seq_lens) < 3:
            print(f'Skipping sequence {s}, not enough peaks')
            continue
        
        Apeaks = sorted(peaks_dic_list[0][s], key=lambda p: p.start)
        
        for Bpeaks in peaks_dic_list[1:]:
            
            common_peaks[s] = []
            
            Bpeaks = sorted(Bpeaks[s], key=lambda p: p.start)
            # Merge two peaks touching each other for each dataset
            Bpeaks = merge_ranges(Bpeaks)
            Apeaks = merge_ranges(Apeaks)
            
            a, b = next(Apeaks), next(Bpeaks)

            while True:

                # Find wether the two current peaks intersect
                max_start = max(a.start, b.start)
                min_stop = min(a.stop, b.stop)
                
                # Try to increment the range with the earlier 
                # stopping value:
                try:
                    if max_start < min_stop:
                        merged = range(min(a.start, b.start),
                                       max(a.stop, b.stop))
                        if a.stop < b.stop:
                            b = merged
                            a = next(Apeaks)
                        else:
                            a = merged
                            b = next(Bpeaks)

                    elif a.stop <= b.start:
                        common_peaks[s].append(a)
                        a = next(Apeaks)

                    else:
                        common_peaks[s].append(b)
                        b = next(Bpeaks)

                except StopIteration:
                    break
            
            Apeaks = list(common_peaks[s])
    
    total_n = sum([len(common_peaks[s]) for s in common_peaks])
    print(f'Found a total of {total_n} common peaks')
    return common_peaks



def rebuild_atac_adatas(adatas, common, layer=None):
    
    with tqdm(total=len(adatas)*len(common), 
              position=0, leave=True) as pbar:
        for i, ann in enumerate(adatas):
            
            peaks = generate_peakdic(ann.var.index.tolist())
            merg = {}
            print(f'Anndata shape: {ann.shape}')
            present = []
            obs_save = ann.obs
            
            start = time.time()
            
            for seq in common:
                
                print(seq)
                
                peaks[seq] = sorted(peaks[seq], key=lambda p: p.start)
                common[seq] = sorted(common[seq], key=lambda p: p.start)
                #print(f'{seq} - {len(list(common[seq]))}')

                # Merge two peaks touching each other for each dataset
                peaks[seq] = merge_ranges(peaks[seq])
                commonset = merge_ranges(common[seq])
            
                a, c = next(peaks[seq]), next(commonset)

                while True:
                    
                    nexta, nextc = False, False
                    a_str = range2peak(seq, a)
                    c_str = range2peak(seq, c)

                    # Find wether the two current peaks intersect
                    int_start = max(a.start, c.start)
                    int_stop = min(a.stop, c.stop)
                    
                    # If it's the same peak
                    if c == a:
                        present.append(a_str)
                        nexta, nextc = True, True

                    # If not the same peak but they intersect
                    elif int_start < int_stop:
                        if c_str in merg.keys():
                            merg[c_str].append(a_str)
                        else:
                            merg[c_str] = [a_str]

                        if a.stop < c.stop:
                            nexta = True
                        else:
                            nextc = True

                    # If common peak not found in this anndata, add
                    elif a.stop <= c.start: 
                        nexta = True
                    else:
                        nextc = True
                        
                    try:
                        if nexta:
                            a = next(peaks[seq])
                        if nextc:
                            c = next(commonset)
                            
                    except StopIteration:
                        try:
                            a = next(peaks[seq])
                        except StopIteration:
                            break
                
                pbar.update(1)
                
            print(f'For loop: {time.time() - start} seconds')
            start = time.time()
            
            feat_types = ann.var['feature_types'][0]
            genome = ann.var['genome'][0]

            # Change index name for merged peaks with 1 value
            old2new = {merg.pop(p)[0]:p for p in merg.copy() 
                       if len(merg[p]) == 1}
            old_merged = list(old2new.keys())
            
            # Otherwise create new anndata for merged peaks with
            # more than 1 old peak
            new_merg = list(merg.keys())

            Xmerg = [ann[:,merg[p]].X.sum(axis=1) for p in merg]
            merg_ann = ad.AnnData(
                X=np.asarray(Xmerg).squeeze().T, 
                obs=ann.obs, 
                var=pd.DataFrame({
                    'gene_ids': new_merg, 
                    'feature_types':feat_types, 
                    'genome':genome}, 
                    index=new_merg
                )
            )
            
            print(f'Creating merge matrix: {time.time() - start} seconds')
            start = time.time()
            
            # Copy the peaks that are the same in the commons and
            # the one to merge with only 1 value and translate
            tokeep = present + old_merged
            tokeep = list(set(ann.var_names).intersection(tokeep))
            ann = ann[:,tokeep]
            ann.var.rename(index=old2new, inplace=True)

            # Concatenate own peaks, not present peaks and merged peaks
            ann = ad.concat([ann, merg_ann], axis=1)
            ann.var_names_make_unique()
            
            print(f'Finalizing merge concatenation: {time.time() - start} seconds')
            start = time.time()

            this_peaks = ann.var.index.tolist()
            common_list = [range2peak(s, p) for s in common 
                                      for p in common[s]]
            
            not_present = list(set(common_list) - set(this_peaks))
            not_present_ann = ad.AnnData(
                X=np.zeros([ann.shape[0], len(not_present)]),
                obs=ann.obs, 
                var=pd.DataFrame({
                    'gene_ids': not_present, 
                    'feature_types':feat_types, 
                    'genome':genome}, 
                    index=not_present
                )
            )

            ann = ad.concat([ann, not_present_ann], axis=1)
            ann.var_names_make_unique()
            
            print(f'Adding not present peaks: {time.time() - start} seconds')
            start = time.time()
            
            this_peaks = ann.var.index.tolist()
            final = list(set(this_peaks).intersection(common_list))
            ann = ann[:,final]
            ann.obs = obs_save
            
            print(f'Ending: {time.time() - start} seconds')
            
            print(f'Final anndata shape: {ann[:,final].shape}')
            adatas[i] = ann
