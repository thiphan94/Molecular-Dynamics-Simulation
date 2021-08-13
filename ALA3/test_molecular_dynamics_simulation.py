import unittest
import numpy as np
from molecular_dynamics_simulation import count_frame, calculate_angle, vector


class Test(unittest.TestCase):
    def test_count_frame(self):
        self.assertEqual(count_frame("G\nG\n1\n2"), 2)

    def test_vector(self):
        coordinates = [
            ["1.796", "1.545", "1.229"],
            ["1.758", "1.522", "1.370"],
            ["1.608", "1.509", "1.388"],
            ["1.557", "1.398", "1.442"],
        ]
        vectors = vector(coordinates)
        np.testing.assert_array_almost_equal(
            vectors,
            np.array(
                [
                    [-0.038, -0.023, 0.141],
                    [0.15, 0.013, -0.018],
                    [-0.051, -0.111, 0.054],
                ]
            ),
        )

    def test_calculate_angle(self):
        input_vectors = (
            np.array([-0.038, -0.023, 0.141]),
            np.array([-0.15, -0.013, 0.018]),
            np.array([-0.051, -0.111, 0.054]),
        )
        expected = 121.88490434193339
        np.testing.assert_almost_equal(
            calculate_angle(input_vectors[0], input_vectors[1], input_vectors[2]),
            expected,
        )


if __name__ == "__main__":
    unittest.main()
