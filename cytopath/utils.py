import anndata
import numpy as np
import pandas as pd
from scipy import sparse, stats

from collections import OrderedDict
from typing import Union

import numba

from .markov_chain_sampler import sampling

# Cosine similarity function
@numba.jit(nopython=True, fastmath=True)
def cosine_similarity_numba(u: np.array, v: np.array):
    assert(u.shape[0] == v.shape[0])

    uv = 0
    uu = 0
    vv = 0
    for j in range(u.shape[0]):
        uv += u[j]*v[j]
        uu += u[j]*u[j]
        vv += v[j]*v[j]
    cos_theta = 1
    if uu!=0 and vv!=0:
        cos_theta = uv/np.sqrt(uu*vv)

    return cos_theta

@numba.jit(nopython=True, fastmath=True)
def cosine_similarity(u: np.array, v: np.array):
    #assert(u.shape[0] == v.shape[0])
    out = np.empty(v.shape[0])
    
    for i in range(out.shape[0]):
        out[i] = cosine_similarity_numba(u, v[i,:])
        
    return out

# Function to generate terminal state frequency using undirected sinulations
def undirected_simulations(data, matrix_key = 'T_forward', max_steps=10, sim_number=None,
                           normalize=False, unique=True, num_cores=1, copy=False):

    """Markov sampling of cell sequences starting from all cells based on a cell-cell transition probability matrix.
    Arguments
    ---------
    adata: :class:`~anndata.AnnData`
        Annotated data matrix with transition probabilities, end points and clustering.
    matrix_key: `string`  (default: 'T_forward')
        Key for accessing transition matrix stored in adata.uns
    max_steps: `integer` (default: 10)
        Maximum number steps of simulation.
    sim_number: `integer`(default: None)
        Number of simulations to generate. If less than number of cells, then a subsampling of cells is used as root states. By default equals number of cells.
    normalize: 'Boolean' (default: False)
        Toggle row sum normalization.
    unique: 'Boolean' (default:True)
        Only keep unique samples.
    num_cores: 'integer' (default:1)
        Number of cpu cores to use.
    copy: 'Boolean' (default: False)
        Create a copy of the anndata object or work inplace.
    
    Returns
    -------
    adata.uns["undirected_samples"]["cell_sequences"]: List of arrays containing the index of the cells in a sequence.
    adata.uns["undirected_samples"]["transition_score"]: List with all of the transition probabilities.
    adata.uns["undirected_samples"]["cluster_sequences"]: The sequence of clusters to which the cells belong.
    adata.obs['log_terminal_state_freq']: Frequency of being terminal state of a simulation for each cell.
    """
    
    adata = data.copy()
    
    adata.obs['undirected_sim_included'] = '1'
    adata.obs['undirected_sim_included'] = adata.obs['undirected_sim_included'].astype('category')

    if sim_number is None:
        sim_number=adata.shape[0]

    if sim_number < adata.shape[0]:
        root_cells = np.random.choice(np.arange(adata.shape[0]), size=sim_number, replace=False)
    else:
        root_cells = np.arange(adata.shape[0])

    sampling(adata, matrix_key=matrix_key, normalize=normalize, unique=unique, 
             root_cells=root_cells, end_points=np.arange(adata.shape[0]), 
             cluster_key='undirected_sim_included', auto_adjust=False, min_clusters=1,
             sim_number=sim_number, traj_number=sim_number, max_steps=max_steps, 
             num_cores=num_cores)

    cells, counts = np.unique(adata.uns['samples']['cell_sequences'][:, -1], return_counts=True)

    adata.obs['log_terminal_state_freq'] = float(np.nan)
    adata.obs.iloc[cells, adata.obs.columns.get_loc('log_terminal_state_freq')] = np.log1p(counts/sim_number)

    data = data.copy() if copy else data
    data.obs['log_terminal_state_freq'] = adata.obs['log_terminal_state_freq']
    data.uns['undirected_samples'] = adata.uns['samples']

    if copy: return data

def differentiation_potential(data, scale=True, copy=False):

    adata = adata.copy() if copy else data
    
    adata.obs['differentiation_potential'] = np.apply_along_axis(stats.entropy, 1, adata.uns['trajectories']['cell_fate_probability'].T)
    if scale: adata.obs['differentiation_potential'] = adata.obs['differentiation_potential']/adata.obs['differentiation_potential'].max()

    if copy: return adata

# https://github.com/thomasmaxwellnorman/perturbseq_demo/blob/master/perturbseq/cell_cycle.py
# https://dynamo-release.readthedocs.io/en/v0.95.2/_modules/dynamo/preprocessing/cell_cycle.html

def norm(x, **kwargs):
    """calculate the norm of an array or matrix"""
    if sparse.issparse(x):
        return sp.linalg.norm(x, **kwargs)
    else:
        return np.linalg.norm(x, **kwargs)

def einsum_correlation(X, Y_i, type="pearson"):
    """calculate pearson or cosine correlation between X (genes/pcs/embeddings x cells) and the velocity vectors Y_i
    for gene i"""

    if type == "pearson":
        X -= X.mean(axis=1)[:, None]
        Y_i -= np.nanmean(Y_i)
    elif type == "cosine":
        X, Y_i = X, Y_i
    elif type == "spearman":
        # check this
        X = stats.rankdata(X, axis=1)
        Y_i = stats.rankdata(Y_i)
    elif type == "kendalltau":
        corr = np.array([stats.kendalltau(x, Y_i)[0] for x in X])
        return corr[None, :]

    X_norm, Y_norm = norm(X, axis=1), norm(Y_i)

    if Y_norm == 0:
        corr = np.zeros(X_norm.shape[0])
    else:
        corr = np.einsum("ij, j", X, Y_i) / (X_norm * Y_norm)[None, :]

    return corr

def group_corr(adata: anndata.AnnData, layer: Union[str, None], gene_list: list) -> tuple:
    """Measures the correlation of all genes within a list to the average expression of all genes within that
    list (used for cell cycle position calling)

    Arguments
    ---------
        adata: :class:`~anndata.AnnData`
            an anndata object.
        layer: `str` or None
            The layer of data to use for calculating correlation. If None, use adata.X.
        gene_list: list of gene names

    Returns
    ---------
        (valid_gene_list, corr): A tuple of valid gene names and the correlation coefficient of each gene with
        the mean expression of all.
    """

    # returns list of correlations of each gene within a list of genes with the total expression of the group
    tmp = adata.var_names.intersection(gene_list)
    # get the location of gene names
    intersect_genes = [adata.var.index.get_loc(i) for i in tmp]

    if len(intersect_genes) == 0:
        raise Exception(f"your adata doesn't have any gene from the gene_list {gene_list}.")

    if layer is None:
        expression_matrix = adata.X[:, intersect_genes]
    else:
        expression_matrix = adata.layers[layer][:, intersect_genes]

    avg_exp = expression_matrix.mean(axis=1)
    cor = (
        einsum_correlation(
            np.array(expression_matrix.A.T, dtype="float"),
            np.array(avg_exp.A1, dtype="float"),
        )
        if sparse.issparse(expression_matrix)
        else einsum_correlation(
            np.array(expression_matrix.T, dtype="float"),
            np.array(avg_exp, dtype="float"),
        )
    )

    # get back to gene names again
    return np.array(adata.var.index[intersect_genes]), cor.flatten()

def refine_gene_list(
    adata: anndata.AnnData,
    layer: Union[str, None],
    gene_list: list,
    threshold: Union[float, None],
    return_corrs: bool = False,
) -> list:
    """Refines a list of genes by removing those that don't correlate well with the average expression of
    those genes

    Parameters
    ----------
        adata: :class:`~anndata.AnnData`
            an anndata object.
        layer: `str` or None
            The layer of data to use for calculating correlation. If None, use adata.X.
        gene_list: list of gene names
        threshold: threshold on correlation coefficient used to discard genes (expression of each gene is
            compared to the bulk expression of the group and any gene with a correlation coefficient less
            than this is discarded)
        return_corrs: whether to return the correlations along with the gene names (default: False)

    Returns
    -------
        Refined list of genes that are well correlated with the average expression trend
    """

    gene_list, corrs = group_corr(adata, layer, gene_list)
    if return_corrs:
        return corrs[corrs >= threshold]
    else:
        return gene_list[corrs >= threshold]


def group_score(adata: anndata.AnnData, layer: Union[str, None], gene_list: list):
    """Scores cells within population for expression of a set of genes. Raw expression data are first
    log transformed, then the values are summed, and then scores are Z-normalized across all cells.

    Arguments
    ---------
        adata: :class:`~anndata.AnnData`
            an anndata object.
        layer: `str` or None
            The layer of data to use for calculating correlation. If None, use adata.X.
        gene_list: list of gene names

    Returns
    -------
        Z-scored expression data
    """

    tmp = adata.var_names.intersection(gene_list)
    # use indices
    intersect_genes = [adata.var_names.get_loc(i) for i in tmp]

    if len(intersect_genes) == 0:
        raise Exception(f"your adata doesn't have any gene from the gene_list {gene_list}.")

    if layer is None:
        expression_matrix = adata.X[:, intersect_genes]
    else:
        expression_matrix = adata.layers[layer][:, intersect_genes]

    # TODO FutureWarning: Index.is_all_dates is deprecated, will be removed in a future version.
    # check index.inferred_type instead
    if layer is None or layer.startswith("X_"):
        scores = expression_matrix.sum(1).A1 if sparse.issparse(expression_matrix) else expression_matrix.sum(1)
    else:
        if sparse.issparse(expression_matrix):
            expression_matrix.data = np.log1p(expression_matrix.data)
            scores = expression_matrix.sum(1).A1
        else:
            scores = np.log1p(expression_matrix).sum(1)

    scores = (scores - scores.mean()) / scores.std()

    return scores


def batch_group_score(adata: anndata.AnnData, layer: Union[str, None], gene_lists: list) -> OrderedDict:
    """Scores cells within population for expression of sets of genes. Raw expression data are first
    log transformed, then the values are summed, and then scores are Z-normalized across all cells.
    Returns an OrderedDict of each score.

    Arguments
    ---------
        adata: an anndata object.
        layer: `str` or None
            The layer of data to use for calculating correlation. If None, use adata.X.
        gene_lists: list of lists of gene names

    Returns
    -------
        an OrderedDict of each score.
    """

    batch_scores = OrderedDict()
    for gene_list in gene_lists:
        batch_scores[gene_list] = group_score(adata, layer, gene_lists[gene_list])
    return batch_scores


def get_cell_phase_genes(
    adata: anndata.AnnData,
    layer: Union[str, None],
    refine: bool = True,
    threshold: Union[float, None] = 0.3,
) -> list:
    """Returns a list of cell-cycle-regulated marker genes, filtered for coherence

    Arguments
    ---------
        adata: an anndata object.
        layer: `str` or None (default: `None`)
            The layer of data to use for calculating correlation. If None, use adata.X.
        refine: `bool` (default: `True`)
            whether to refine the gene lists based on how consistent the expression is among
            the groups
        threshold: `float` or None (default: `0.3`)
            threshold on correlation coefficient used to discard genes (expression of each
            gene is compared to the bulk expression of the group and any gene with a correlation
            coefficient less than this is discarded)

    Returns
    -------
        a list of cell-cycle-regulated marker genes that show strong co-expression
    """

    cell_phase_genes = OrderedDict()
    cell_phase_genes["G1-S"] = pd.Series(
        [
            "ARGLU1",
            "BRD7",
            "CDC6",
            "CLSPN",
            "ESD",
            "GINS2",
            "GMNN",
            "LUC7L3",
            "MCM5",
            "MCM6",
            "NASP",
            "PCNA",
            "PNN",
            "SLBP",
            "SRSF7",
            "SSR3",
            "ZRANB2",
        ]
    )
    cell_phase_genes["S"] = pd.Series(
        [
            "ASF1B",
            "CALM2",
            "CDC45",
            "CDCA5",
            "CENPM",
            "DHFR",
            "EZH2",
            "FEN1",
            "HIST1H2AC",
            "HIST1H4C",
            "NEAT1",
            "PKMYT1",
            "PRIM1",
            "RFC2",
            "RPA2",
            "RRM2",
            "RSRC2",
            "SRSF5",
            "SVIP",
            "TOP2A",
            "TYMS",
            "UBE2T",
            "ZWINT",
        ]
    )
    cell_phase_genes["G2-M"] = pd.Series(
        [
            "AURKB",
            "BUB3",
            "CCNA2",
            "CCNF",
            "CDCA2",
            "CDCA3",
            "CDCA8",
            "CDK1",
            "CKAP2",
            "DCAF7",
            "HMGB2",
            "HN1",
            "KIF5B",
            "KIF20B",
            "KIF22",
            "KIF23",
            "KIFC1",
            "KPNA2",
            "LBR",
            "MAD2L1",
            "MALAT1",
            "MND1",
            "NDC80",
            "NUCKS1",
            "NUSAP1",
            "PIF1",
            "PSMD11",
            "PSRC1",
            "SMC4",
            "TIMP1",
            "TMEM99",
            "TOP2A",
            "TUBB",
            "TUBB4B",
            "VPS25",
        ]
    )
    cell_phase_genes["M"] = pd.Series(
        [
            "ANP32B",
            "ANP32E",
            "ARL6IP1",
            "AURKA",
            "BIRC5",
            "BUB1",
            "CCNA2",
            "CCNB2",
            "CDC20",
            "CDC27",
            "CDC42EP1",
            "CDCA3",
            "CENPA",
            "CENPE",
            "CENPF",
            "CKAP2",
            "CKAP5",
            "CKS1B",
            "CKS2",
            "DEPDC1",
            "DLGAP5",
            "DNAJA1",
            "DNAJB1",
            "GRK6",
            "GTSE1",
            "HMG20B",
            "HMGB3",
            "HMMR",
            "HN1",
            "HSPA8",
            "KIF2C",
            "KIF5B",
            "KIF20B",
            "LBR",
            "MKI67",
            "MZT1",
            "NUF2",
            "NUSAP1",
            "PBK",
            "PLK1",
            "PRR11",
            "PSMG3",
            "PWP1",
            "RAD51C",
            "RBM8A",
            "RNF126",
            "RNPS1",
            "RRP1",
            "SFPQ",
            "SGOL2",
            "SMARCB1",
            "SRSF3",
            "TACC3",
            "THRAP3",
            "TPX2",
            "TUBB4B",
            "UBE2D3",
            "USP16",
            "WIBG",
            "YWHAH",
            "ZNF207",
        ]
    )
    cell_phase_genes["M-G1"] = pd.Series(
        [
            "AMD1",
            "ANP32E",
            "CBX3",
            "CDC42",
            "CNIH4",
            "CWC15",
            "DKC1",
            "DNAJB6",
            "DYNLL1",
            "EIF4E",
            "FXR1",
            "GRPEL1",
            "GSPT1",
            "HMG20B",
            "HSPA8",
            "ILF2",
            "KIF5B",
            "KPNB1",
            "LARP1",
            "LYAR",
            "MORF4L2",
            "MRPL19",
            "MRPS2",
            "MRPS18B",
            "NUCKS1",
            "PRC1",
            "PTMS",
            "PTTG1",
            "RAN",
            "RHEB",
            "RPL13A",
            "SRSF3",
            "SYNCRIP",
            "TAF9",
            "TMEM138",
            "TOP1",
            "TROAP",
            "UBE2D3",
            "ZNF593",
        ]
    )

    if refine:
        for phase in cell_phase_genes:
            cur_cell_phase_genes = None
            if adata.var_names[0].isupper():
                cur_cell_phase_genes = cell_phase_genes[phase]
            elif adata.var_names[0][0].isupper() and adata.var_names[0][1:].islower():
                cur_cell_phase_genes = [gene.capitalize() for gene in cell_phase_genes[phase]]
            else:
                cur_cell_phase_genes = [gene.lower() for gene in cell_phase_genes[phase]]

            cell_phase_genes[phase] = refine_gene_list(adata, layer, cur_cell_phase_genes, threshold)

    return cell_phase_genes


def get_cell_phase(
    adata: anndata.AnnData,
    layer: str = None,
    gene_list: Union[OrderedDict, None] = None,
    refine: bool = True,
    threshold: Union[float, None] = 0.3,
) -> pd.DataFrame:
    """Compute cell cycle phase scores for cells in the population

    Arguments
    ---------
        adata: :class:`~anndata.AnnData`
        layer: `str` or None (default: `None`)
            The layer of data to use for calculating correlation. If None, use adata.X.
        gene_list: `OrderedDict` or None (default: `None`)
            OrderedDict of marker genes to use for cell cycle phases. If None, the default
            list will be used.
        refine: `bool` (default: `True`)
            whether to refine the gene lists based on how consistent the expression is among
            the groups
        threshold: `float` or None (default: `0.3`)
            threshold on correlation coefficient used to discard genes (expression of each
            gene is compared to the bulk expression of the group and any gene with a correlation
            coefficient less than this is discarded)

    Returns
    -------
        Cell cycle scores indicating the likelihood a given cell is in a given cell cycle phase
    """

    # get list of genes if one is not provided
    if gene_list is None:
        cell_phase_genes = get_cell_phase_genes(adata, layer, refine=refine, threshold=threshold)
    else:
        cell_phase_genes = gene_list

    adata.uns["cell_phase_genes"] = cell_phase_genes
    # score each cell cycle phase and Z-normalize
    phase_scores = pd.DataFrame(batch_group_score(adata, layer, cell_phase_genes))
    normalized_phase_scores = phase_scores.sub(phase_scores.mean(axis=1), axis=0).div(phase_scores.std(axis=1), axis=0)

    normalized_phase_scores_corr = normalized_phase_scores.transpose()
    normalized_phase_scores_corr["G1-S"] = [1, 0, 0, 0, 0]
    normalized_phase_scores_corr["S"] = [0, 1, 0, 0, 0]
    normalized_phase_scores_corr["G2-M"] = [0, 0, 1, 0, 0]
    normalized_phase_scores_corr["M"] = [0, 0, 0, 1, 0]
    normalized_phase_scores_corr["M-G1"] = [0, 0, 0, 0, 1]

    phase_list = ["G1-S", "S", "G2-M", "M", "M-G1"]

    # final scores for each phaase are correlation of expression profile with vectors defined above
    cell_cycle_scores = normalized_phase_scores_corr.corr()
    tmp = -len(phase_list)
    cell_cycle_scores = cell_cycle_scores[tmp:].transpose()[: -len(phase_list)]

    # pick maximal score as the phase for that cell
    cell_cycle_scores["cell_cycle_phase"] = cell_cycle_scores.idxmax(axis=1)
    cell_cycle_scores["cell_cycle_phase"] = cell_cycle_scores["cell_cycle_phase"].astype("category")
    cell_cycle_scores["cell_cycle_phase"].cat.set_categories(phase_list, inplace=True)

    def progress_ratio(x, phase_list):
        ind = phase_list.index(x["cell_cycle_phase"])
        return x[phase_list[(ind - 1) % len(phase_list)]] - x[phase_list[(ind + 1) % len(phase_list)]]

    # interpolate position within given cell cycle phase
    cell_cycle_scores["cell_cycle_progress"] = cell_cycle_scores.apply(
        lambda x: progress_ratio(x, list(phase_list)), axis=1
    )
    cell_cycle_scores.sort_values(
        ["cell_cycle_phase", "cell_cycle_progress"],
        ascending=[True, False],
        inplace=True,
    )

    # order of cell within cell cycle phase
    cell_cycle_scores["cell_cycle_order"] = cell_cycle_scores.groupby("cell_cycle_phase").cumcount()
    cell_cycle_scores["cell_cycle_order"] = cell_cycle_scores.groupby("cell_cycle_phase")["cell_cycle_order"].apply(
        lambda x: x / (len(x) - 1)
    )

    return cell_cycle_scores

def cell_cycle_scores(
    adata: anndata.AnnData,
    layer: Union[str, None] = None,
    gene_list: Union[OrderedDict, None] = None,
    refine: bool = True,
    threshold: float = 0.3,
    copy: bool = False,
) -> anndata.AnnData:
    """Call cell cycle positions for cells within the population. If more direct control is desired,
    use get_cell_phase.

    Arguments
    ---------
        adata: an anndata object.
        layer: `str` or None (default: `None`)
            The layer of data to use for calculating correlation. If None, use adata.X.
        gene_list: OrderedDict of marker genes to use for cell cycle phases. If None, the default
            list will be used.
        refine: `bool` (default: `True`)
            whether to refine the gene lists based on how consistent the expression is among
            the groups
        threshold: `float` or None (default: `0.3`)
            threshold on correlation coefficient used to discard genes (expression of each
            gene is compared to the bulk expression of the group and any gene with a correlation
            coefficient less than this is discarded)
        copy:
            If true, copy the original AnnData object and return it

    Returns
    -------
        Returns an updated adata object with cell_cycle_phase as new column in .obs and a new data
        frame with `cell_cycle_scores` key to .obsm where the cell cycle scores indicating the likelihood a
        given cell is in a given cell cycle phase.
    """
    adata = adata.copy() if copy else adata


    cell_cycle_scores = get_cell_phase(
        adata,
        layer=layer,
        refine=refine,
        gene_list=gene_list,
        threshold=threshold,
    )

    cell_cycle_scores.index = adata.obs_names[cell_cycle_scores.index.values.astype("int")]
    adata.obs["cell_cycle_phase_prediction"] = cell_cycle_scores["cell_cycle_phase"].astype("category")
    adata.obsm["cell_cycle_scores"] = cell_cycle_scores.loc[adata.obs_names, :]
