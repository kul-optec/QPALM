import os

from conan import ConanFile
from conan.errors import ConanInvalidConfiguration
from conan.tools.cmake import CMakeToolchain, CMake, cmake_layout
from conan.tools.build import can_run


class QPALMRecipe(ConanFile):
    name = "qpalm"
    version = "1.2.6"
    package_type = "library"

    # Optional metadata
    license = "LGPL-3.0-or-later"
    author = "Pieter P <pieter.p.dev@outlook.com>"
    url = "https://github.com/kul-optec/QPALM"
    description = "Proximal Augmented Lagrangian method for Quadratic Programs"
    topics = ("optimization", "quadratic-program", "qp", "alm")

    # Binary configuration
    settings = "os", "compiler", "build_type", "arch"
    bool_qpalm_options = {
        "with_cxx": True,
        "with_python": False,
        "with_julia": False,
        "with_fortran": False,
        "with_mex": False,
    }
    options = (
        {
            "shared": [True, False],
            "fPIC": [True, False],
        }
        | {k: [True, False] for k in bool_qpalm_options}
        | {"with_conan_python": [True, False]}
    )
    default_options = (
        {
            "shared": False,
            "fPIC": True,
        }
        | bool_qpalm_options
        | {"with_conan_python": False}
    )

    # Sources are located in the same place as this recipe, copy them to the recipe
    exports_sources = (
        "CMakeLists.txt",
        "QPALM/*",
        "LADEL/*",
        "examples/*",
        "test/*",
        "LICENSE",
        "README.md",
    )

    generators = ("CMakeDeps",)

    def requirements(self):
        self.requires("ladel/0.0.4", transitive_headers=True)
        self.test_requires("gtest/1.17.0")
        if self.options.with_cxx:
            self.requires("eigen/3.4.0", transitive_headers=True)
        if self.options.with_python:
            self.requires("pybind11/3.0.1")
            if self.options.with_conan_python:
                self.requires("tttapa-python-dev/3.13.7")

    def build_requirements(self):
        self.tool_requires("cmake/[>=4.1 <5]")

    def configure(self):
        if self.options.get_safe("with_mex"):
            self.options["ladel/*"].with_mex = True

    def config_options(self):
        if self.settings.get_safe("os") == "Windows":
            self.options.rm_safe("fPIC")

    def validate(self):
        if self.options.with_python and not self.options.with_cxx:
            msg = "Python interface requires C++. Set 'with_cxx=True'."
            raise ConanInvalidConfiguration(msg)

    def layout(self):
        cmake_layout(self)

    def generate(self):
        tc = CMakeToolchain(self)
        tc.variables["QPALM_WITH_EXAMPLES"] = False
        for k in self.bool_qpalm_options:
            value = getattr(self.options, k, None)
            if value is not None and value.value is not None:
                tc.variables["QPALM_" + k.upper()] = bool(value)
        if can_run(self):
            tc.variables["QPALM_FORCE_TEST_DISCOVERY"] = True
        tc.generate()

    def build(self):
        cmake = CMake(self)
        cmake.configure()
        cmake.build()
        cmake.test()

    def package(self):
        cmake = CMake(self)
        cmake.install()

    def package_info(self):
        self.cpp_info.set_property("cmake_find_mode", "none")
        self.cpp_info.set_property("cmake_file_name", "QPALM")
        self.cpp_info.builddirs.append(os.path.join("lib", "cmake", "QPALM"))
        if self.options.with_cxx:
            self.cpp_info.builddirs.append(os.path.join("lib", "cmake", "QPALM_cxx"))
        if self.options.with_fortran:
            self.cpp_info.builddirs.append(
                os.path.join("lib", "cmake", "QPALM_fortran")
            )
