#!/usr/bin/env python

import blade_coverage
import argparse
import db
from os import environ
import sys
from os.path import join
from turbine_config import TurbineConfig
from turbine import Turbine
import xml.etree.ElementTree as ET

def jusante_grids():
    """ Return the int numbers from jusante grids.

    Returns:
        grid numbers

    Examples:
        >>> grid_nums = jusante_grids()
    """

    grid_nums = range(0, 15)
    grid_nums.extend(range(17, 20))
    grid_nums.extend(range(22, 25))
    grid_nums.extend(range(67, 70))
    grid_nums.extend(range(72, 77))
    grid_nums.extend(range(78, 80))
    return grid_nums

def montante_grids():
    """ Return the int numbers from montante grids.

    Returns:
        grid numbers

    Examples:
        >>> grid_nums = montante_grids()
    """

    grid_nums = range(30, 50)
    grid_nums.extend(range(51, 55))
    grid_nums.extend(range(56, 60))
    grid_nums.extend([77])
    return grid_nums


def lip_grids():
    """ Return the int numbers from lip grids.

    Returns:
        grid numbers

    Examples:
        >>> grid_nums = lip_grids()
    """

    return [0, 1, 2]


def border_grids():
    """ Return the int numbers from border grids.

    Returns:
        grid numbers

    Examples:
        >>> grid_nums = border_grids()
    """

    grid_nums = range(60, 65)
    grid_nums.extend(range(25, 29))
    grid_nums.extend([85])
    return grid_nums


area_db = {'jusante': join(environ['PYTHON_DATABASE'],'FACE'),
           'montante': join(environ['PYTHON_DATABASE'],'FACE'),
           'lip': join(environ['PYTHON_DATABASE'],'LIP'),
           'border': join(environ['PYTHON_DATABASE'],'BORDER')
          }

grids_available = {'jusante': jusante_grids,
                   'montante': montante_grids,
                   'lip': lip_grids,
                   'border': border_grids
                  }


def coverage(args):
    turbine = True
    DB = db.DB(area_db[args.Area], turbine)
    psa, _ = DB.get_rail_configuration_n(args.Grids, args.Config_file, 'sum', args.ans)
    print psa
    return

def grids_out(args):
    print grids_available[args.Area]()
    return

def areas_out(args):
    print area_db.keys()
    return
    
def validate(args):
    config = TurbineConfig.load(args.Config_file)
    config.environment.load = args.Environment_file # Set given path as environment.xml path

    env_file = ET.parse(config.environment.load)
    robot_xml = env_file.find('Robot')
    robot_xml.attrib['file'] = args.Robot_file
    env_file.write(config.environment.load) # Set given path as robot.xml path

    turbine = Turbine(config)
    DB = db.DB(join(environ['PYTHON_DATABASE'],area_db[args.Area]),turbine)
    success = blade_coverage.base_grid_validation_parser(turbine, DB, args.Grid, args.Trajectory_file)
    if not success:
        print "could not find the required coating trajectory"
        sys.exit(1)

parser = argparse.ArgumentParser(description='Coating script.')

subparsers = parser.add_subparsers(help='Functions')

# Function to return PSAlpha
parser_p = subparsers.add_parser('position', help='Compute the rail and robot base position.')
parser_p.add_argument('Config_file', type=str, help='Config file path - turbine.cfg.')
parser_p.add_argument('Area', choices=area_db.keys(), type=str, help='Area to be coated.')
parser_p.add_argument('Grids', nargs='+', type=int, help='Grids to be coated.')
parser_p.add_argument('-ans', type=int, help='Answer number.', default=0)
parser_p.set_defaults(func=coverage)

# Function to return Areas
parser_a = subparsers.add_parser('areas', help='Return available areas.')
parser_a.set_defaults(func=areas_out)


# Function to return Grids
parser_g = subparsers.add_parser('grids', help='Return available grids given area.')
parser_g.add_argument('Area', choices=area_db.keys(), type=str, help='Area to be coated.')
parser_g.set_defaults(func=grids_out)

# Function to return Trajectory File
parser_j = subparsers.add_parser('plan', help='Return the trajectory file for grid coating.')
parser_j.add_argument('Area', choices=area_db.keys(), type=str, help='Area to be coated.')
parser_j.add_argument('Grid', type=int, help='Grid to be coated.')
parser_j.add_argument('Config_file', type=str, help='Config file path - turbine.cfg')
parser_j.add_argument('Environment_file', type=str, help='Environment file path - env.xml')
parser_j.add_argument('Robot_file', type=str, help='Robot file path - robot.xml')
parser_j.add_argument('Trajectory_file', type=str, help='Path to save traj.xml')
parser_j.set_defaults(func=validate)

args = parser.parse_args()
args.func(args)
