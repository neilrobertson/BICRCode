import scipy as sp
import scipy.stats as stats


def pvalue(lmeth,lunmeth,rmeth,runmeth):
    fisherInputTable = [[lmeth, lunmeth], [rmeth, runmeth]]
    oddsRatio, pvalue = stats.fisher_exact(fisherInputTable, "two-sided")
    return pvalue