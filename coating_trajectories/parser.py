import blade_coverage
import argparse

area_db = {'jusante':['db','jusante'],
           'montante': ['db','montante'],
           'lip': ['db_lip','lip']
          }

grids_available = {'jusante':blade_coverage.jusante_grids,
                   'montante': blade_coverage.montante_grids,
                   'lip': blade_coverage.lip_grids
                  }

def coverage(args):
    db_directories, grid_path = area_db[args.Area]
    print blade_coverage.robot_base_position(db_directories, grid_path, args.Grid)

def grids_out(args):
    print grids_available[args.Area]()

parser = argparse.ArgumentParser(description='Coating script.')
parser.add_argument('Area', choices=area_db.keys(), type=str,
                    help = 'Area to be coated.')

subparsers = parser.add_subparsers(help='Functions')

# Function to return PSAlpha
parser_p = subparsers.add_parser('p', help='Compute the rail and robot base position.')
parser_p.add_argument('Grid', type=int, help='Grid to be coated.')
parser_p.set_defaults(func=coverage)

# Function to return Grids
parser_g = subparsers.add_parser('g', help='Return available grids given area.')
parser_g.set_defaults(func=grids_out)


args = parser.parse_args()
args.func(args)
