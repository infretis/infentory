import numpy as np
from infretis.classes.orderparameter import OrderParameter
import MDAnalysis as mda
from MDAnalysis.analysis.distances import distance_array
from MDAnalysis.transformations.wrap import make_whole
from MDAnalysis.analysis import rms
import torch

class NeuralNetDiffmap(OrderParameter):
    def __init__(self):
        super().__init__(description="An order parameter for chignolin.")
        # indices of the protein
        self.protein = [i for i in range(166)]
        # create a universe, used for the rmsd calculation
        self.u = mda.Universe("average-protein-structure.pdb")
        # store reference structure positions
        self.ref = self.u.select_atoms('backbone').atoms.positions
        # the neural network takes in a set of distances between
        # protein backbone atoms, the indices are contained in the
        # file distances
        tmp = np.loadtxt("distances.txt")
        self.dist = tmp[:,:2].astype(int)
        self.features = (self.dist[:,0],self.dist[:,1])
        self.mean = tmp[:,2]
        self.std = tmp[:,3]
        # load torch model
        self.X = torch.zeros(1,self.dist.shape[0])
        self.model = torch.load("model-diffusion-map.pt")
        self.model.eval()

    def calculate(self, system):
        # retrieve the positions
        self.u.atoms.positions = system.pos[self.protein]*10
        # calculate the distance array
        dist_arr = distance_array(self.u.atoms.positions, 
                             self.u.atoms.positions,
                             box=self.u.dimensions)
        # make whole for rmsd calculation
        make_whole(self.u.atoms)
        rmsd = rms.rmsd(self.u.select_atoms("backbone").positions,
                self.ref,superposition=True)
        # pick out the features to give to the neural network
        # subtract the mean and divide by the std. dev
        # (these where calculated during training)
        x=(dist_arr[(self.features)]-self.mean)/self.std
        self.X[0]=torch.from_numpy(x)
        # evaluate the model
        H=self.model(self.X)[0]
        H1,H2=H[0].item(),H[1].item()
        # find the value s where the point (H1,H2) on the diffusion map 
        # is equal to the the point on a parabola given by
        # (a*x^2+s,x). The solution is s = H1 - a*H2^2.
        # (Set H1=a*x^2+s and H2 = x. => H1 = a*H2^2+s => S = H1 - a*H2^2
        a=1/15
        diffmap_order = H1-a*H2**2

        # write some additional info
        return [diffmap_order, H1, H2, rmsd]
