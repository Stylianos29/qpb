# QCD Package for Brillouin improved fermions

## Dependencies
To build and use this project, you **need** the following dependencies
installed:
- **OpenMPI** (`mpicc`)
- **C-lime** ([GitHub link](https://github.com/usqcd-software/c-lime))
- **GNU Autotools**
- **GSL** (GNU Scientific Library)
- **SymPy** Python library

---

## Setup Instructions

### Step 1: Clone QPB Locally
1. Clone this repository to a directory of your choice, e.g., `qpb_build/`.
   ```bash
   git clone <qpb-repository-link> qpb_build/
   ```
2. Navigate to a directory of your choice **other than and outside**
   `qpb_build/qpb/` (the root directory of the QPB project) for the next steps.

### Step 2: Install C-lime
1. Clone the C-lime repository:
   ```bash
   git clone https://github.com/usqcd-software/c-lime.git
   cd c-lime
   ```
2. Build and install C-lime: 
* Run:
   ```bash
   ./autogen.sh
   ```
   * If successful, this command will generate additional files, including a
     `configure` script.
   ### Troubleshooting:
   * If you see an error like `autoreconf: command not found`, install GNU
     Autotools.
   * If you encounter `macro ... is obsolete`, run:
   ```bash
   autoupdate
   ```
   * Once the issue is resolved, rerun:
   ```bash
   ./autogen.sh
   ```
* Choose an installation directory within `qpb/`, such as a new directory like
  `qpb_build/qpb/install/`, and set the `LIMEPREFIX` variable. **NOTE:** Avoid
  using `qpb_build/qpb/` as the installation directory to prevent conflicts and
  the risk of overwriting important files. It’s recommended to build in a
  separate directory from the source to maintain a clean and organized project
  structure.
   ```bash
   export LIMEPREFIX=<full-path-to-install-directory>
   ```
* Configure and build using the MPI compiler `mpicc` (do not use `gcc` to avoid
  compatibility issues with the `qpb` project):
   ```bash
   ./configure --prefix=$LIMEPREFIX CC=mpicc CFLAGS=-O3
   make
   make install
   ```

### Step 3: Configure QPB
1. Navigate to the QPB directory:
   ```bash
   cd qpb_build/qpb
   ```
2. Choose a `Makefile.in.<...>` file (e.g., `Makefile.in.CYCLONE`) and link it:
   ```bash
   ln -s Makefile.in.CYCLONE Makefile.in
   ```
   **Note:** Ensure `CFLAGS` and `LDFLAGS` in this file point correspondingly to
   the **full paths** of:
   * `qpb/install/include/`
   * `qpb/install/lib/`

### Step 4: Build QPB
1. Build the *library*:
   ```bash
   cd qpb/lib
   make
   ```
   #### Troubleshooting:
   If you see `fatal error: gsl/gsl_rng.h: No such file or directory`, install
   GSL.
2. Build *main programs*:
   ```bash
   cd qpb/mainprogs
   make
   ```  
   #### Troubleshooting:
   Errors or warnings during either step are unacceptable and require restarting
   the setup process from the beginning. Before restarting:
   1. **Verify Dependencies**: Ensure all required dependencies are installed
      and properly linked to the Makefiles.
   2. **Timestamp Conflicts**: If issues persist, updating timestamps before
      rebuilding the library may help:
      ```bash
      cd qpb/lib
      touch *
      ```  
      **Note:** This step has been observed to resolve issues potentially
      related to newer SymPy versions, but the exact cause may vary.
   3. Once these checks are complete, repeat the setup process from the
      beginning.

### Step 5: Running QPB Programs

1. Navigate to a `mainprogs` directory (e.g.,
   `overlap-kl/ginsparg-wilson-relation`) and modify the `params.ini` file as
   needed. For more information on filling in the parameters file, consult the
   guide (N/A yet).  
   <!-- TODO: Guide for filling in the parameters files -->
2. Run the program with MPI (if supported). For example, using `mpirun`:
   ```bash
   mpirun -n 8 --bind-to core --report-bindings ./ginsparg-wilson-relation geom=2,2,2 params.ini
   ```
   Alternatively, you can use `srun`:
   ```bash
   srun -n 8 ./ginsparg-wilson-relation geom=2,2,2 params.ini
   ```
   