
from collections import defaultdict

class EpitopeCollection(object):
    def __init__(self, binding_predictions):
        self.binding_predictions = list(sorted(set(binding_predictions)))

    def __len__(self):
        return len(self.binding_predictions)

    def __iter__(self):
        return iter(self.binding_predictions)

    def __str__(self):
        return "<EpitopeCollection with %d elements>" % (len(self),)

    def __repr__(self):
        return str(self)

    def filter(self, filter_fn):
        return self.__class__(x for x in self if filter_fn(x))

    def filter_by_binding_threshold(self, threshold):
        return self.filter(lambda x: x.measure.is_binder(x.value, threshold))

    def filter_by_percentile_rank(self, max_rank):
        return self.filter(lambda x: x.percentile_rank <= max_rank)

    def groupby(self, key_fn):
        groups = defaultdict(list)
        for binding_prediction in self.binding_predictions:
            key = key_fn(binding_prediction)
            groups[key].append(binding_prediction)
        # want to create an EpitopeCollection for each group
        # but need to write the construction in terms of
        # self.__class__ so that this works with derived classes
        return {
            key: self.__class__(binding_predictions)
            for (key, binding_predictions)
            in groups.items()
        }

    def groupby_allele(self):
        return self.groupby(key_fn=lambda x: x.allele)

    def groupby_peptide(self):
        return self.groupby(key_fn=lambda x: x.peptide)

    def groupby_allele_and_peptide(self):
        return self.groupby(key_fn=lambda x: (x.allele, x.peptide))
