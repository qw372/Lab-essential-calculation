import numpy as np
import matplotlib.pyplot as plt

from dataclasses import dataclass

@dataclass
class Ray:
    """
    Ray

    origin: np.ndarray, origin of ray
    direction: np.ndarray, direction of ray
    """

    origin: np.ndarray
    direction: np.ndarray

    def __post_init__(self) -> None:
        """
        Post-initialization, check vector lengths
        """
        assert len(self.origin) == 3
        assert len(self.direction) == 3

    def norm(self) -> None:
        """
        Normalize the direction vector
        """
        self.direction /= np.linalg.norm(self.direction)

class PointLightSource:
    """
    Point Light source
    """

    def __init__(self, origin: np.ndarray, direction: np.ndarray) -> None:
        """
        origin: np.ndarray, origin of light source
        """
        
        self.origin = origin
        self.direction = direction

    def generate_rays(self):
        """
        Generate rays from the light source

        num_rays: number of rays to generate
        """

        return Ray(self.origin, self.direction)
    
