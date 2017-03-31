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

parser = argparse.ArgumentParser(description='Coating script.')
group = parser.add_mutually_exclusive_group()

# Functions
group.add_argument('-p', '--position', nargs=2, metavar=('Area', 'Grid'),
                    help='Compute the rail and robot base position. Area = '+str(area_db.keys())+'. Grid = int.')
group.add_argument('-g', '--grids', nargs=1, metavar=('Area'), 
                    help='Return available grids given area . Area = '+str(area_db.keys())+'.')


args = parser.parse_args()

if args.position:
    db_directories, grid_path = area_db[args.position[0]]
    print blade_coverage.robot_base_position(db_directories, grid_path, int(args.position[1]))
elif args.grids:
    print grids_available[args.grids[0]]()
