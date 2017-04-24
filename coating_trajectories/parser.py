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
           'lip': 'LIP'
          }

grids_available = {'jusante':blade_coverage.jusante_grids,
                   'montante': blade_coverage.montante_grids,
                   'lip': blade_coverage.lip_grids
                  }

def coverage(args):
    DB = db.DB(area_db[args.Area],turb)
    psa = DB.get_rail_configuration_n(args.Grids, 'sum', 0)
    print psa

def grids_out(args):
    print grids_available[args.Area]()

parser = argparse.ArgumentParser(description='Coating script.')

subparsers = parser.add_subparsers(help='Functions')

# Function to return PSAlpha
parser_p = subparsers.add_parser('position', help='Compute the rail and robot base position.')
parser_p.add_argument('Area', choices=area_db.keys(), type=str, help='Area to be coated.')
parser_p.add_argument('Grids', nargs='+', type=int, help='Grids to be coated.')
parser_p.set_defaults(func=coverage)

# Function to return Grids
parser_g = subparsers.add_parser('grids', help='Return available grids given area.')
parser_g.add_argument('Area', choices=area_db.keys(), type=str, help='Area to be coated.')
parser_g.set_defaults(func=grids_out)


args = parser.parse_args()
args.func(args)
