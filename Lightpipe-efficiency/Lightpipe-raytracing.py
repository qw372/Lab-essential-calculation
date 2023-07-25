import numpy as np
import matplotlib.pyplot as plt

from dataclasses import dataclass

def Gram_Schmidt_ortho(a: np.ndarray, b: np.ndarray, b_normalized: bool = False) -> np.ndarray:
    """
    return a's component that's orthogonal to b
    """

    assert a.shape == b.shape

    a = a.astype(float)
    b = b.astype(float)

    if not b_normalized:
        b /= np.linalg.norm(b)

    return a - np.dot(a, b) * b


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

        assert self.origin.shape == (3,)
        assert self.direction.shape == (3,)

        self.origin = self.origin.astype(float)
        self.direction = self.direction.astype(float)

    def norm(self) -> None:
        """
        Normalize the direction vector
        """

        self.direction /= np.linalg.norm(self.direction)


class AbstractLightSource:
    """
    AbstractLightSource class, used as base class for other light sources
    """
    
    def __init__(self):
        pass

    def generate_rays(self, number_of_rays: int) -> list[Ray]:
        pass


class PointLightSource(AbstractLightSource):
    """
    Point light source
    """

    def __init__(self, origin: np.ndarray, axis: np.ndarray, field_of_view: float) -> None:
        """
        Define a ligth source

        origin: np.ndarray, origin of light source
        axis: np.ndarray, symmetry axis of light source
        field_of_view: float, field of view of light source (full angle)
        """
        
        self.origin = origin
        self.axis = axis
        self.field_of_view = field_of_view

        assert self.origin.shape == (3,)
        assert self.axis.shape == (3,)
        assert self.field_of_view > 0 and self.field_of_view <= 2 * np.pi

        self.origin = self.origin.astype(float)
        self.axis = self.axis.astype(float)
        self.field_of_view = float(self.field_of_view)

        self.axis /= np.linalg.norm(self.axis) # normalize vector

    def generate_rays(self, number_of_rays: int) -> list[Ray]:
        """
        Generate rays from the light source

        number_of_rays: number of rays to generate
        """

        theta_range = self.field_of_view / 2
        rng = np.random.default_rng(123)

        # theta has probability density function sin(theta)
        # cos(theta) has uniform distribution
        cos_theta_random = rng.uniform(np.cos(theta_range), 1, size=number_of_rays)
        theta_random = np.arccos(cos_theta_random)

        # phi has uniform distribution
        phi_random = rng.uniform(0, 2 * np.pi, size=number_of_rays)

        # generate a unit vector has the same direction as the light source axis's shortest component
        # so it won't be parallel to the light source axis
        vec = np.array([0, 0, 0])
        vec[np.argmin(np.abs(self.axis))] = 1

        # generate two unit vectors that are orthogonal to the light source axis, using cross product
        unit_1 = np.cross(self.axis, vec)
        unit_1 /= np.linalg.norm(unit_1)
        unit_2 = np.cross(unit_1, self.axis)

        # ray directions are calculated from its projection on the light source axis and the two orthogonal unit vectors
        rays = [Ray(origin=self.origin, direction=self.axis*np.cos(theta)+unit_1*np.sin(theta*np.cos(phi))+unit_2*np.sin(theta)*np.sin(phi)) for theta, phi in zip(theta_random, phi_random)]
        
        return rays


class AbstractOpticComponent:
    """
    AbstractOpticComponent class, used as base class for other optic components
    """
    
    def __init__(self):
        pass

    def intersect(self, ray: Ray):
        pass

    def plot(self, ax: plt.axes):        
        pass

class AbsorbingPlane(AbstractOpticComponent):
    """
    Defines a (infinitely large) plane that absorbs all light incident on it
    """
    
    def __init__(self, origin: np.ndarray, normal: np.ndarray) -> None:
        """ 
        origin: np.ndarray, origin of plane
        normal: np.ndarray, normal vector of plane
        """

        self.origin = origin
        self.normal = normal

        assert self.origin.shape == (3,)
        assert self.normal.shape == (3,)

        self.origin = self.origin.astype(float)
        self.normal = self.normal.astype(float)

    def intersect(self, ray: Ray) -> tuple:
        """
        Find the hit point of a ray on the plane

        ray: Ray, incident ray
        """

        if np.dot(ray.direction, self.normal) < 1e-8:
            # ray is (quasi-)parallel to plane
            a = None
            hit_point = None

        else:
            # the end point of the ray can be parametrized as ray.origin + a * ray.direction
            # the end point of the ray must lie on the plane, i.e. np.dot(ray.origin + a * ray.direction - self.origin, self.normal) = 0
            a = np.dot(self.origin - ray.origin, self.normal) / np.dot(ray.direction, self.normal)

            # only cosider forward going rays
            if a > 0:
                hit_point = ray.origin + a * ray.direction
            else:
                a = None
                hit_point = None

        # refelcted/refracted ray direction, it's None in this case
        ray_dir_next = None

        return (a, hit_point, ray_dir_next)
    
    def plot(self, ax: plt.axes) -> None:
        """
        Plot the optic component

        ax: plt.axes, axes to plot on
        """

        pass

    
class ReflectingCylindricalSurface(AbstractOpticComponent):
    """
    Defines (infinitely long) cylindrical surface that refelcs all the light incident on it
    """

    def __init__(self, origin: np.ndarray, axis: np.ndarray, radius: float) -> None:
        """
        origin: np.ndarray, origin of cylinder
        axis: np.ndarray, axis of cylinder
        radius: float, radius of cylinder
        """

        self.origin = origin
        self.axis = axis
        self.radius = radius

        assert self.origin.shape == (3,)
        assert self.axis.shape == (3,)

        self.origin = self.origin.astype(float)
        self.axis = self.axis.astype(float)
        self.radius = float(self.radius)
        assert self.radius > 0

        self.axis /= np.linalg.norm(self.axis) # normalize vector

    def intersect(self, ray: Ray) -> tuple:
        """
        Find the hit point and refelcted ray direction of a ray incident on the cylinder surface

        ray: Ray, incident ray
        """
        
        if np.linalg.norm(Gram_Schmidt_ortho(ray.direction, self.axis, b_normalized=True)) < 1e-8:
            # gram-schmidt orthogonalization to find the component of the ray direction that is orthogonal to the cylinder axis
            # if it's zero, ray is (quasi-)parallel to cylinder axis
            a = None
            hit_point = None
            reflected_ray_dir = None

        else:
            # Ignore the component of the ray direction that is parallel to the cylinder axis
            # reduce all directional vectors to the component that is orthogonal to the cylinder axis
            # reduce the problem to a 2D plane that's orthogonal to the cylinder axis

            # vector from cylinder origin to ray origin
            ray_origin_to_cyl = self.origin - ray.origin
            ray_origin_to_cyl_2D = Gram_Schmidt_ortho(ray_origin_to_cyl, self.axis, b_normalized=True)
            ray_origin_to_cyl_2D_norm = np.linalg.norm(ray_origin_to_cyl_2D)

            # reduced vector of ray direction
            ray_dir_2D = Gram_Schmidt_ortho(ray.direction, self.axis, b_normalized=True)

            # cos of the angle between 1)ray direction and 2)ray origin to cylinder origin
            # on this 2D plane, ray origin, cylinder origin and hit point form a triangle
            # use law of cosines to find the distance between ray origin and hit point
            if ray_origin_to_cyl_2D_norm == 0:
                # ray origin is on the cylinder axis
                cos_angle = 0
            else:
                cos_angle = np.dot(ray_dir_2D, ray_origin_to_cyl_2D) / (np.linalg.norm(ray_dir_2D) * ray_origin_to_cyl_2D_norm)

            # alpha is the discriminant of the quadratic equation
            alpha = self.radius**2 - ray_origin_to_cyl_2D_norm**2 + (ray_origin_to_cyl_2D_norm * cos_angle)**2

            if alpha <= 0:
                # ray misses cylinder
                a = None
                hit_point = None
                reflected_ray_dir = None

            else:
                # ray hits cylinder
                beta_minus = ray_origin_to_cyl_2D_norm * cos_angle - np.sqrt(alpha)
                beta_plus = ray_origin_to_cyl_2D_norm * cos_angle + np.sqrt(alpha)

                if beta_minus >= 0:
                    # ray hits cylinder from outside
                    a = beta_minus / np.linalg.norm(ray_dir_2D)
                    hit_point = ray.origin + a * ray.direction # in 3D

                else:
                    # ray hits cylinder from inside
                    a = beta_plus / np.linalg.norm(ray_dir_2D)
                    hit_point = ray.origin + a * ray.direction

                # norm to the cylinder surface at the hit point in 2D plane
                cyl_surface_norm_2D = ray_origin_to_cyl_2D - a * ray_dir_2D

                # component of ray direction that's perpendicular to the cylinder surface, in 2D
                ray_dir_2D_perpendicular = Gram_Schmidt_ortho(ray_dir_2D, cyl_surface_norm_2D, b_normalized=False)

                # component of ray direction that's parallel to the cylinder surface, in 2D
                ray_dir_2D_parallel = ray_dir_2D - ray_dir_2D_perpendicular

                # reflected ray direction in 2D, ray_dir_2D's component that's perpendicular to cylinder axis is preserved, the parallel component is flipped
                reflected_ray_dir_2D = ray_dir_2D_perpendicular + ray_dir_2D_parallel * (-1)

                # reflected ray direction in 3D, add back the component that's parallel to cylinder axis
                reflected_ray_dir = reflected_ray_dir_2D + (ray.direction - ray_dir_2D)

        return (a, hit_point, reflected_ray_dir)
    
    def plot(self, ax: plt.axes) -> None:
        """
        Plot the optic component

        ax: plt.axes, axes to plot on
        """

        pass

    
class Container(AbstractOpticComponent):
    """
    Defines a container that contains other optic components. 
    Use as absorbing boundary condition for all rays so they don't propagate forever
    """
    
    def __init__(self):
        pass

    def intersect(self, ray: Ray):
        pass

    def plot(self, ax: plt.axes):
        pass

class OpticalSystem:
    """
    An Optical system includes light source, optic components and containers, as well as functions to trace rays through the system
    """
    
    def __init__(self):
        self.light_source = []
        self.optical_components = []
        self.containers = None

        self.fig = plt.figure()
        self.ax = self.fig.add_subplot(projection='3d')
        self.ax.set_xlabel("x")
        self.ax.set_ylabel("y")
        self.ax.set_zlabel("z")

    def add_light_source(self, light_source) -> None:
        """
        Add a light source to the system
        """
        
        self.light_source.append(light_source)

    def add_optic_component(self, optic_component: AbstractOpticComponent) -> None:
        """
        Add an optic component to the system
        """

        self.optical_components.append(optic_component)

    def add_container(self, container: Container) -> None:
        """
        Add a container to the system
        """
        
        self.add_container = container

    def generate_rays(self, number_of_rays: np.ndarray) -> None:
        """ 
        Generate rays from all light sources

        number_of_rays: np.ndarray, number of rays to generate for each light source
        """
        
        assert len(self.light_source) > 0
        assert len(self.light_source) == len(number_of_rays)

        self.rays = []
        for light_source, number_of_rays in zip(self.light_source, number_of_rays):
            self.rays.extend(light_source.generate_rays(number_of_rays))

    def propagate_rays(self, max_num_refelction: float = 50) -> None:
        """
        Trace rays

        max_num_refelction: float, maximum number of refelctions allowed for each ray
        """
        
        cloest_hit_point_list_total = []
        for ray in self.rays:
            ray_dir_next = np.array([0, 0, 0])
            cloest_hit_point_list = ray.origin
            i = 0

            while i < max_num_refelction:
                # find the closest hit point and the corresponding optic component
                a_list = []
                hit_point_list = []
                ray_dir_next_list = []

                for optic_component in self.optical_components:
                    a, hit_point, ray_dir_next = optic_component.intersect(ray)

                    if a is not None:
                        a_list.append(a)
                        hit_point_list.append(hit_point)
                        ray_dir_next_list.append(ray_dir_next)

                if len(a_list) == 0:
                    # no hit point found, ray does't intersect with any optic component
                    break

                else:
                    # find the closest hit point
                    a_min_index = np.argmin(a_list)
                    hit_point_min = hit_point_list[a_min_index]
                    ray_dir_next_min = ray_dir_next_list[a_min_index]

                    # update ray origin and direction
                    ray.origin = hit_point_min
                    ray.direction = ray_dir_next_min

                    # record hit point
                    cloest_hit_point_list = np.vstack((cloest_hit_point_list, ray.origin))

                if ray.direction is None:
                    break

                # propagate ray along its direction by a small distance to avoid self-intersection
                ray.origin += ray.direction * 1e-6

            cloest_hit_point_list_total.append(cloest_hit_point_list)

        self.hit_point_list = cloest_hit_point_list_total

    def plot_trace(self) -> None:
        """
        Plot the rays and optical components
        """

        for ray_trace in self.hit_point_list:
            self.ax.plot(ray_trace[:, 0], ray_trace[:, 1], ray_trace[:, 2], "r-")

    def plot_optic_compoents(self) -> None:
        """
        Plot the optical components
        """

        for optic_component in self.optical_components:
            optic_component.plot(self.ax)

    def plot_container(self) -> None:
        """
        Plot the containers
        """

        self.container.plot(self.ax)

    def show_plot(self) -> None:
        """
        Show the plot
        """

        plt.show()


# define light source
src = PointLightSource(origin=np.array([0, 0, 0]), axis=np.array([1, 1, 0]), field_of_view=np.pi/6)

# define optic components
cyl = ReflectingCylindricalSurface(origin=np.array([0, 0, 0]), axis=np.array([1, 0, 0]), radius=5)
pln = AbsorbingPlane(origin=np.array([10, 0, 0]), normal=np.array([1, 0, 0]))

# define optical system
sys = OpticalSystem()
sys.add_light_source(src)
sys.add_optic_component(cyl)
sys.add_optic_component(pln)

# trace rays
sys.generate_rays(number_of_rays=np.array([100]))
sys.propagate_rays(max_num_refelction=10)

# plot
sys.plot_trace()
sys.plot_optic_compoents()
sys.show_plot()

# To-do:
# 1. plot optical components
# 2. implement container
# 3. implement more light sources and optic components