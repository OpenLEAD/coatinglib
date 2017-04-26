import blade_coverage
import argparse
import db
from os import environ
from os.path import realpath, join
from turbine_config import TurbineConfig
from turbine import Turbine

dir_test = join(environ['PYTHON_COAT'],'test')
environ['OPENRAVE_DATA'] = str(dir_test)
cfg = TurbineConfig.load('turbine_unittest.cfg',dir_test)
turb = Turbine(cfg)

area_db = {'jusante':'FACE',
           'montante': 'FACE',
           'lip': 'LIP',
           'border':'BORDER'
          }

grids_available = {'jusante':blade_coverage.jusante_grids,
                   'montante': blade_coverage.montante_grids,
                   'lip': blade_coverage.lip_grids,
                   'border': blade_coverage.border_grids
                  }


def coverage(args):
    DB = db.DB(area_db[args.Area],turb)
    psa, fcomp = DB.get_rail_configuration_n(args.Grids, 'sum', args.ans)
    print psa
    return

def grids_out(args):
    print grids_available[args.Area]()
    return

def validate(args):
    DB = db.DB(area_db[args.Area],turb)
    score, joints = blade_coverage.base_grid_validation(DB, args.Grid)
    print score
    return 

parser = argparse.ArgumentParser(description='Coating script.')

subparsers = parser.add_subparsers(help='Functions')

# Function to return PSAlpha
parser_p = subparsers.add_parser('position', help='Compute the rail and robot base position.')
parser_p.add_argument('Area', choices=area_db.keys(), type=str, help='Area to be coated.')
parser_p.add_argument('Grids', nargs='+', type=int, help='Grids to be coated.')
parser_p.add_argument('-ans', type=int, help='Answer number.', default=0)
parser_p.set_defaults(func=coverage)

# Function to return Grids
parser_g = subparsers.add_parser('grids', help='Return available grids given area.')
parser_g.add_argument('Area', choices=area_db.keys(), type=str, help='Area to be coated.')
parser_g.set_defaults(func=grids_out)

# Function to return Grids
parser_j = subparsers.add_parser('joints', help='Return joints for grid coating.')
parser_j.add_argument('Area', choices=area_db.keys(), type=str, help='Area to be coated.')
parser_j.add_argument('Grid', type=int, help='Grid to be coated.')
parser_j.set_defaults(func=validate)

args = parser.parse_args()
args.func(args)
