import unittest
import turbine

class TestTurbine(unittest.TestCase):

    @classmethod
    def setUpClass(cls):
        super(TestTurbine, cls).setUpClass()
        cls.turb = turbine.Turbine("turbine_unittest.cfg",False)

    def test_load_wrong(self):
        with self.assertRaises(turbine.ConfigFileError):
            turbwrong = turbine.Turbine("wrong",False)
        
    # environment
    def test_parsed_load(self):
        self.assertEqual(TestTurbine.turb.environment.load,"../Turbina/env_mh12_0_16.xml")

    def test_parsed_z_floor_level(self):
        self.assertEqual(TestTurbine.turb.environment.z_floor_level,-3.22)

    def test_parsed_primary_safe_margin(self):
        self.assertEqual(TestTurbine.turb.environment.primary_safe_margin,0.2)

    def test_parsed_secondary_safe_margin(self):
        self.assertIs(type(TestTurbine.turb.environment.secondary_safe_margin), float)
        self.assertEqual(TestTurbine.turb.environment.secondary_safe_margin,2)

    def test_parsed_robot_level_difference(self):
        self.assertEqual(TestTurbine.turb.environment.robot_level_difference,0.01)

    def test_parsed_blade_angle(self):
        self.assertIs(type(TestTurbine.turb.environment.blade_angle),float)
        self.assertEqual(TestTurbine.turb.environment.blade_angle,0)

    def test_parsed_rotor_angle(self):
        self.assertEqual(TestTurbine.turb.environment.rotor_angle,-5.73)

    def test_parsed_x_max(self):
        self.assertEqual(TestTurbine.turb.environment.x_max,1.5)

    def test_parsed_x_min(self):
        self.assertEqual(TestTurbine.turb.environment.x_min,-1)

    def test_parsed_y_max(self):
        self.assertEqual(TestTurbine.turb.environment.y_max,2)

    def test_parsed_y_min(self):
        self.assertEqual(TestTurbine.turb.environment.y_min,-2)

    def test_parsed_rail_angle_limit(self):
        self.assertEqual(TestTurbine.turb.environment.rail_angle_limit,0.34906585)

    def test_parsed_rail_angle_limit(self):
        self.assertEqual(TestTurbine.turb.environment.rail_angle_mean,0.523598776)

    # coating
    def test_parsed_min_distance(self):
        self.assertEqual(TestTurbine.turb.coating.min_distance,0.21)

    def test_parsed_ideal_distance(self):
        self.assertEqual(TestTurbine.turb.coating.ideal_distance,0.23)

    def test_parsed_max_distance(self):
        self.assertEqual(TestTurbine.turb.coating.max_distance,0.24)

    def test_parsed_angle_tolerance(self):
        self.assertEqual(TestTurbine.turb.coating.angle_tolerance,30)

    def test_parsed_coating_speed(self):
        self.assertEqual(TestTurbine.turb.coating.coating_speed,0.667)

    def test_parsed_parallel_gap(self):
        self.assertEqual(TestTurbine.turb.coating.parallel_gap, 0.003)
        
    # model
    def test_parsed_nose_radius(self):
        self.assertEqual(TestTurbine.turb.model.nose_radius,0.59)

    def test_parsed_runner_radius(self):
        self.assertEqual(TestTurbine.turb.model.runner_radius,3.75)

if __name__ == '__main__':
    unittest.main()
