from sklearn.metrics import f1_score
import sparseconvnet as scn
import torch.utils.data
from tqdm import tqdm
import numpy as np
import importlib
import warnings
import argparse
import torch
import glob
import math
import sys
import os

warnings.filterwarnings('ignore')
BASE_DIR = os.path.dirname(os.path.abspath(__file__))
ROOT_DIR = BASE_DIR
CLASS_LABELS = ['0', '1']
model_path = 'model'

aa_list = ['A', 'R', 'N', 'D', 'C', 'Q', 'E', 'G', 'H', 'I', 'L', 'K', 'M', 'F', 'P', 'S', 'T', 'W', 'Y', 'V', 'X']
atom_list = ['N', 'C', 'O', 'S', 'H', 'X']

aa_classes = {}
for idx, aa in enumerate(aa_list):
    aa_classes[aa] = idx
atom_classes = {}
for idx, atom in enumerate(atom_list):
    atom_classes[atom] = idx

label2color = {'0': [0, 255, 255], '1': [255, 0, 0]}

def one_hot(length, position):
    zeros = [0 for _ in range(length)]
    zeros[position] = 1
    return zeros


def parse_args():
    '''PARAMETERS'''
    parser = argparse.ArgumentParser('Model')
    parser.add_argument('--gpu', type=str, default='0', help='GPU to use [default: GPU 0]')
    parser.add_argument('--output', type=str, required=True,help='Output direction [e.g /data/output/]')
    parser.add_argument('--data', type=str, required=True,help='Folder include lig, atom [e.g /data/COACH420_data/]')
    parser.add_argument('--select_list', type=str, required=True,help='Path for selected list [e.g /data/list]')
    parser.add_argument('--num_vote', type=int, default=25, help='Number of voting [default: 25]')
    parser.add_argument('--seed', type=int, default=1, help='Random Seed for voting [default: 1]')

    return parser.parse_args()


def process_pc_xyz(protein):
    '''
    Only can preocess .xyz data

    columns of new_data :
    coords: xyz
    feature:
        0-20: aa
        21-26: atom
        27: seq_idx
    y: lebel
    raw_coord: raw coordinate
    seq: residue
    '''

    protein_data = []
    seq = []
    with open(protein) as f:
        for line in f:
            line = line.strip().split()
            seq.append(line[0])
            aa_onehot = one_hot(21, aa_classes[line[1][0]])
            atom_onehot = one_hot(6, atom_classes[line[1][1]])
            xyz = [float(line[axis]) for axis in range(2, 5)]
            point_feature = xyz + aa_onehot + atom_onehot
            seg_label = [line[5]]
            line_data = point_feature + seg_label
            protein_data.append(line_data)
    protein_data = np.array(protein_data)
    'Normalize'
    coords = protein_data[:, 0:3].astype(float)
    raw_coord = coords
    feature = protein_data[:, 3:-1].astype(int)
    y = protein_data[:, -1].astype(float)
    centroid = np.mean(coords, axis=0)
    coords = coords - centroid
    m = np.max(np.sqrt(np.sum(coords ** 2, axis=1)))
    coords = coords / m

    return (coords, feature, y, raw_coord, seq)


def main(args):
    sys.path.append(model_path)
    DATA_PATH = args.data
    val_file = DATA_PATH + '/*_atom.xyz'
    SELECT_LIST = args.select_list
    val_file = glob.glob(val_file)

    def predict_binding_site(chain_coords, chain_feature, random_seed=1):
        vote_pool = torch.zeros(chain_coords.shape[0], 2)
        vote_num = torch.zeros(chain_coords.shape[0], 2)

        chain_coords_news = []
        chain_features = []
        chain_point_id = []

        '''Make Inference'''
        for se in range(args.num_vote):
            seed = random_seed * se
            m = np.eye(3)
            m[0][0] *= np.random.randint(0, 2) * 2 - 1
            m *= scale
            np.random.seed(seed)
            theta = np.random.rand() * 2 * math.pi
            m = np.matmul(m, [[math.cos(theta), math.sin(theta), 0], [-math.sin(theta), math.cos(theta), 0],
                              [0, 0, 1]])
            np.random.seed(seed)
            chain_coords_new = np.matmul(chain_coords, m) + full_scale / 2 + np.random.uniform(-2, 2, 3)
            Min = chain_coords_new.min(0)
            Max = chain_coords_new.max(0)
            np.random.seed(seed)
            offset = - Min + np.clip(full_scale - Max + Min - 0.001, 0, None) * np.random.rand(3) + np.clip(
                full_scale - Max + Min + 0.001, None, 0) * np.random.rand(3)
            chain_coords_new += offset
            idxs = (chain_coords_new.min(1) >= 0) * (chain_coords_new.max(1) < full_scale)
            coords = torch.Tensor(chain_coords_new[idxs]).long()
            chain_coords_news.append(torch.cat([coords, torch.LongTensor(coords.shape[0], 1).fill_(se)], 1))
            chain_features.append(torch.Tensor(chain_feature[idxs]))
            chain_point_id.append(torch.Tensor(np.nonzero(idxs)[0]))

        chain_coords_news = torch.cat(chain_coords_news, 0)
        chain_features = torch.cat(chain_features, 0).float()
        chain_point_ids = torch.cat(chain_point_id, 0).long()

        if use_cuda:
            chain_features = chain_features.cuda()

        predictions = classifier([chain_coords_news, chain_features])
        predictions = torch.nn.functional.softmax(predictions)
        vote_pool.index_add_(0, chain_point_ids, predictions.cpu())
        vote_num.index_add_(0, chain_point_ids, torch.ones_like(predictions.cpu()))
        vote_pool = vote_pool / vote_num

        return vote_pool

    '''HYPER PARAMETER'''
    use_cuda = True if torch.cuda.is_available() else False
    FEATURE_DIMENTION = 27
    output_dir = args.output
    if not os.path.exists(output_dir):
        os.mkdir(output_dir)

    print('Inference examples: %d' % len(val_file))

    '''MODEL LOADING'''
    MODEL = importlib.import_module('unet')
    num_class = 2
    scale = MODEL.scale
    full_scale = MODEL.full_scale
    classifier = MODEL.Model(FEATURE_DIMENTION, num_class)
    if use_cuda:
        classifier = classifier.cuda()

    print('#classifer parameters: %d' % (sum([x.nelement() for x in classifier.parameters()])))
    try:
        print('Load model %s' % model_path)
        if use_cuda:
            classifier.load_state_dict(torch.load(os.path.join(model_path,'scale_80.pth')))
        else:
            classifier.load_state_dict(torch.load(os.path.join(model_path, 'scale_80.pth'), map_location='cpu'))
    except:
        raise ValueError('Cannot load pretrain model from %s!!!' % model_path)

    with torch.no_grad():
        classifier.eval()
        total_pred = []
        total_label = []
        select_list = []
        with open(SELECT_LIST) as f:
            lines = f.readlines()
            for line in lines:
                select_list.append(line.strip())
        val_file = [i for i in val_file if i.split('/')[-1][0:i.split('/')[-1].find('_atom')] in select_list]
        for i, protein in tqdm(enumerate(val_file), total=len(val_file)):
            pro = protein.split('/')[-1]
            pro = pro[:pro.find('_atom.xyz')]
            coords, feature, label, raw_coord, residue_idx = process_pc_xyz(protein)

            chain_idx = np.array([i.split('|')[0] for i in residue_idx])
            '''Inference in Chain Level'''
            vote_pool = torch.zeros(coords.shape[0], 2)
            chain_ids = torch.zeros(coords.shape[0], 1)
            for ch_i, chain_id in enumerate(np.unique(chain_idx)):
                chain_vote_pool = predict_binding_site(coords[chain_idx == chain_id], feature[chain_idx == chain_id],
                                                       random_seed=args.seed)
                point_ids = torch.Tensor(np.nonzero(chain_idx == chain_id))[0].long()
                vote_pool.index_add_(0, point_ids, chain_vote_pool)
                chain_ids.index_add_(0, point_ids, torch.ones(len(chain_vote_pool)) * ch_i)

            confident_map = vote_pool.data.numpy()[:, 1]
            atom_choose = np.zeros(vote_pool.shape[0])
            atom_choose[confident_map > 0.5] = 1

            pro_base_dir = os.path.join(output_dir, pro)
            if not os.path.exists(pro_base_dir):
                os.mkdir(pro_base_dir)
            visual_dir = os.path.join(pro_base_dir, 'visual')
            if not os.path.exists(visual_dir):
                os.mkdir(visual_dir)

            label[label != 0] = 1
            data_output = np.loadtxt(protein, dtype=str)
            data_output[:,5] = confident_map
            data_output = data_output.astype(str)
            np.savetxt(os.path.join(pro_base_dir, '%s_output.xyz' % pro),data_output,fmt='%s')

            fout = open(os.path.join(visual_dir, '%s_atom_pred.obj' % pro), 'w')
            fout_gt = open(os.path.join(visual_dir, '%s_atom_gt.obj' % pro), 'w')
            for j in range(raw_coord.shape[0]):
                color = label2color[str(int(atom_choose[j]))]
                color_gt = label2color[str(int(label[j]))]
                fout.write('v %f %f %f %d %d %d\n' % (raw_coord[j, 0], raw_coord[j, 1], raw_coord[j, 2],
                                                      color[0], color[1], color[2]))
                fout_gt.write('v %f %f %f %d %d %d\n' % (raw_coord[j, 0], raw_coord[j, 1], raw_coord[j, 2],
                                                         color_gt[0], color_gt[1], color_gt[2]))

            fout.close()
            fout_gt.close()
            total_pred += list(atom_choose)
            total_label += list(label)

        print('Num of protein %d' % len(val_file))
        f1 = f1_score(total_pred, total_label)
        if f1 > 0:
            print('Atom Level F1_score %.3f' % f1)


if __name__ == '__main__':
    args = parse_args()
    os.environ["CUDA_VISIBLE_DEVICES"] = args.gpu
    main(args)
