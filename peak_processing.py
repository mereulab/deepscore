import re
import gc
import pandas as pd
import numpy as np
import matplotlib.pyplot as plt
from tqdm import tqdm
import anndata as ad



def filter_peaks_width(peakdic, width_range):
    """
    Filter a dictionary of peaks by peak range width.
    - peakdic: Dictionary of peaks per sequence (chromosome)
    - with_range: Tuple with the minimum and maximum peak
        peak width to keep. Other peaks will be discarded
    """

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
    """
    Transform a peak from string to Python range object.
    """

    try:
        chrom, r = peak.split(':')
        r = r.split('-')
        rang = range(int(r[0]), int(r[1])+1)
        return chrom, rang
    except ValueError:
        print(f'Cannot process peak {peak}. Please ensure ',
               'that all peaks have the correct format.\nRun ',
               'ensure_peak_format() on all anndata objects.')



def range2peak(chrom, rang):
    """
    Transform a peak from Python range object to string.
    """
    return f'{chrom}:{rang.start}-{rang.stop-1}'



def generate_peakdic(peak_list):
    """
    Generate a dictionary of peak ranges by sequence
    (chromosome) given a list of peaks in string format.
    """

    dic = {}
    for peak in peak_list:

        chrom, rang = peak2range(peak)

        if chrom in dic:
            dic[chrom].append(rang)
        else:
            dic[chrom] = [rang]

    return dic



def ensure_peak_format(ann, sep=[':', '-']):
    """
    Ensure that the AnnData object have a consistent
    peak format for later comparisons and operations.
    """

    peaks = ann.var.index.tolist()
    non_alphanum = "[^0-9a-zA-Z.]+"

    regex = [re.findall(non_alphanum, p) for p in peaks]
    if all(r==sep for r in regex):
        print('Peak format is consistent')
        return ann
    else:
        print('Peak format is inconsistent. Reformatting.')
        peaks = [p.replace(r[0], sep[0], 1).replace(r[1], sep[1], 1)
                for r, p in zip(regex, peaks)]
        ann.var = ann.var.reindex(index=peaks)
        return ann



def merge_ranges(peak_list):
    """
    Given a list of peaks, find and merge the peaks that are
    inmediately adjacent to another, resulting in a single peak.
    """

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
    """
    Keep only the sequences belonging to chromosomes
    """
    return [x for x in chroms if x.startswith('chr')]



def get_common_peaks(adatas, only_standard_chroms=False):
    """
    Find a common set of peaks between a list of AnnData
    objects.
    """

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


def rebuild_atac_datas(adatas, common, layer=None,
                       feat_types="Peaks", genome="GRCh38"):
    """
    Rebuild an scATAC AnnData object using a common set of
    peak features given.
    - adatas: list of AnnData objects to rebuild.
    - common: list of peaks present in all AnnData objects.
    - layer: TODO.
    - feat_types: Type of the features as specified in 10X output.
    - genome: Reference genome.
    """

    with tqdm(total=len(adatas)*len(common),
              position=0, leave=True) as pbar:
        for i, ann in enumerate(adatas):

            print(f'Anndata shape: {ann.shape}')
            peaks = generate_peakdic(ann.var.index.tolist())
            merg = {}
            present = []
            obs_save = ann.obs

            for seq in common:

                peaks[seq] = sorted(peaks[seq], key=lambda p: p.start)
                common[seq] = sorted(common[seq], key=lambda p: p.start)
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

            # Change index name for merged peaks with 1 value
            old2new = {merg.pop(p)[0]:p for p in merg.copy()
                       if len(merg[p]) == 1}
            old_merged = list(old2new.keys())

            # Otherwise create new anndata for merged peaks with
            # more than 1 old peak
            new_merg = list(merg.keys())

            # This is a little mess but it's the only fast way
            # CHAPUCERO, VER SI SE PUEDE ARREGLAR
            Xmerg = []
            not_present = []
            ann = ann.T
            for p in merg:
                try:
                    idxs = [ann.obs.index.get_loc(k) for k in merg[p]]
                    Xmerg.append(ann.chunk_X(idxs).sum(axis=0))
                except KeyError:
                    print('not found?')
                    not_present.extend(merg[p])
                    ann = ann[~np.in1d(ann.obs_names, merg[p])].copy()
                    new_merg.remove(p)

            ann = ann.T

            # This method is double longer in computational time
            #Xmerg = [ann[:,merg[p]].X.sum(axis=1) for p in merg]

            # If there are merged peaks
            if len(merg) > 0:
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

            # Copy the peaks that are the same in the commons and
            # the one to merge with only 1 value and translate
            tokeep = present + old_merged
            tokeep = list(set(ann.var_names).intersection(tokeep))
            ann = ann[:,tokeep]
            ann.var.rename(index=old2new, inplace=True)

            # Concatenate own peaks, not present peaks and merged peaks
            if len(merg) > 0:
                ann = ad.concat([ann, merg_ann], axis=1)
                ann.var_names_make_unique()

            this_peaks = ann.var.index.tolist()
            common_list = [range2peak(s, p) for s in common
                                      for p in common[s]]

            not_present.extend(list(set(common_list) - set(this_peaks)))
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

            this_peaks = ann.var.index.tolist()
            final = list(set(this_peaks).intersection(common_list))
            ann = ann[:,final]
            ann.obs = obs_save

            print(f'Final anndata shape: {ann[:,final].shape}')
            adatas[i] = ann

            del not_present_ann
            if len(merg) > 0:
                del merg_ann
            gc.collect()
