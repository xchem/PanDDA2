from ..interfaces import *

def get_comparators(datasets: Dict[str, DatasetInterface], filters, debug=False):
    print(f'comparators : before : {len([x for x in datasets])}')

    for filter in filters:
        datasets = filter(datasets)
        if debug:
            print(f'comparators : {type(filter)} : {len([x for x in datasets])}')

    return datasets
