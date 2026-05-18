import numpy as np

from ..interfaces import *

class SparseDMap:

    def __init__(self, data: np.array):
        self.data = data

    @classmethod
    def from_xmap(cls, xmap: CrystallographicGridInterface, dframe: DFrameInterface, debug=False):
        xmap_array = np.array(xmap, copy=False)

        if debug:
            x, y, z = dframe.mask.indicies[0], dframe.mask.indicies[1], dframe.mask.indicies[2]
            mask_range = ((np.min(x), np.max(x)), (np.min(y), np.max(y)), (np.min(z), np.max(z)))
            xmap_array_shape = xmap_array.shape
            print(f'Xmap array shape: {xmap_array_shape} vs mask range: {mask_range}')
        data = xmap_array[dframe.mask.indicies]

        return cls(data)

    ...