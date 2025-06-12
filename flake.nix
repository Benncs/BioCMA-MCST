{
  inputs = {
    nixpkgs.url = "github:NixOS/nixpkgs/nixos";
    flake-utils.url = "github:numtide/flake-utils";
  };

  outputs = { self, nixpkgs, flake-utils }:
    flake-utils.lib.eachDefaultSystem (system:
      let
        pkgs = import nixpkgs {
          inherit system;
          config.allowUnfree = true;
        };

        toolsSrc = [ ./tools ./devutils ];
        appsSrc = toolsSrc ++ [ ./apps ];

        pythonEnv = pkgs.python3.withPackages (ps: with ps; [
          pybind11
        ]);

        llvmPackages = pkgs.llvmPackages_20;

        buildInputs = with pkgs; [
          meson
          ninja
          cmake
          pythonEnv
          pkg-config
          llvmPackages.clang
        ];
      in {
        packages.default = pkgs.stdenv.mkDerivation {
          pname = "name";
          version = "0.1.0";

          src = appsSrc;

          nativeBuildInputs = buildInputs;

          # Uncomment if using Meson
          # mesonFlags = [ "-Dinstall_header=true" ];
          # mesonBuildType = "release";
        };

        devShells.default = pkgs.mkShell {
          buildInputs = buildInputs;
          src = appsSrc;
            pure = true;
          shellHook = ''
            export CC=${llvmPackages.clang}/bin/clang
            export CXX=${llvmPackages.clang}/bin/clang++
            echo "Dev environment ready!"
          '';
        };
      }
    );
}
