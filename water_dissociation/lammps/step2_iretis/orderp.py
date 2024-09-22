import numpy as np
from MDAnalysis.analysis import distances
from infretis.classes.orderparameter import OrderParameter

# Assign global variables:
NMOL = 32
NATOMS = NMOL*3
HIDX = [i for i in range(NATOMS) if i % 3 != 0]
OIDX = list(range(0, NATOMS, 3))


class OrderX(OrderParameter):
    """Daniel watar dissociation op."""
    def __init__(self):
        super().__init__()

    def calculate(self, system):
        # set up system variables
        pos = system.pos
        box = [*system.box,90,90,90]

        # key: OIDX, value: dict of HIDXs and dists
        dinfo = {i: {'hidxs': [], 'dists': []} for i in OIDX}

        # calc distance pairs
        for idx in HIDX:
            dists = distances.distance_array(pos[idx],
                                             pos[OIDX],
                                             box=box)[0]
            o_close = OIDX[np.argmin(dists)]
            oh_dist = dists[np.argmin(dists)]
            dinfo[o_close]['hidxs'].append(idx)
            dinfo[o_close]['dists'].append(oh_dist)

        # check do not have isolated O:
        assert len(dinfo) == NMOL

        # calculate num h per o list
        olist = [len(i['hidxs']) for i in dinfo.values()]

        # if no dissociation:
        if all(i == 2 for i in olist):
            orderp = max([max(i['dists']) for i in dinfo.values()])
            hattach = -1
            oattach_min = -1
            oattach_max = -1
            typep = 0

        # we have 30 2's, h3o+, oh-
        elif sum([i == 2 for i in olist]) == NMOL - 2 and set(olist) == set((1,2,3)):
        # elif sum([i == 2 for i in olist]) == len(olist) - 2:
            o_min = OIDX[np.argmin(olist)]
            o_max = OIDX[np.argmax(olist)]
            dists = distances.distance_array(pos[o_min],
                                             pos[dinfo[o_max]['hidxs']],
                                             box=box)[0]
            orderp = min(dists)
            hattach = dinfo[o_max]['hidxs'][np.argmin(dists)]
            oattach_min = o_min
            oattach_max = o_max
            typep = 1
        else:
            o_mins = np.array(OIDX)[np.array(olist)==1]
            o_maxs = np.array(OIDX)[np.array(olist)==3]
            dists = []
            # for o_min, o_max in zip(o_mins, o_maxs):
            hattach = -1
            for o_min in o_mins:
                for o_max in o_maxs:
                    dists_i = distances.distance_array(pos[o_min],
                                                       pos[dinfo[o_max]['hidxs']],
                                                       box=box)[0]
                    if dists and min(dists_i) < min(dists):
                        hattach = dinfo[o_max]['hidxs'][np.argmin(dists_i)]
                    dists += list(dists_i)
                    oattach_min = o_min
                    oattach_max = o_max
            orderp = min(dists)
            typep = 2

        return [orderp, typep, oattach_min, oattach_max, hattach]
