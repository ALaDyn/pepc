import argparse
import numpy as np

parser = argparse.ArgumentParser()
parser.add_argument("ref", help="specify path to reference data", type=str)
parser.add_argument("source", help="specify paths to compare with ref", type=str)
args = parser.parse_args()

ref = np.genfromtxt(args.ref,unpack=False)
source = np.genfromtxt(args.source,unpack=False)

delta_end = ref[-1] - source[-1]
print np.linalg.norm(delta_end[2:5],np.inf)/np.linalg.norm(ref[-1,2:5],np.inf), np.linalg.norm(delta_end[5:8],np.inf)/np.linalg.norm(ref[-1,5:8],np.inf)


# assert len(ref) == len(source)

# delta =  ref-source

# rel_l2norm_x  = []
# rel_infnorm_x = []
# rel_l2norm_v  = []
# rel_infnorm_v = []

# for i in range(len(delta)):
#     rel_l2norm_x.append( np.linalg.norm(delta[i,2:5],2)     /np.linalg.norm(ref[i,2:5],2))
#     rel_l2norm_v.append( np.linalg.norm(delta[i,2:5],2)     /np.linalg.norm(ref[i,2:5],2))
#     rel_infnorm_x.append(np.linalg.norm(delta[i,2:5],np.inf)/np.linalg.norm(ref[i,2:5],np.inf))    
#     rel_infnorm_v.append(np.linalg.norm(delta[i,5:8],np.inf)/np.linalg.norm(ref[i,5:8],np.inf))

# print rel_infnorm_x[-1], rel_infnorm_v[-1]

 
