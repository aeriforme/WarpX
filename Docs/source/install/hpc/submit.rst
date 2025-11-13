.. _building-submit:

SubMIT (MIT)
============

The `SubMIT cluster <https://submit.mit.edu/>`_ is located at MIT.


Introduction
------------

If you are new to this system, **please see the following resources**:

* `SubMIT user guide <https://submit.mit.edu/submit-users-guide/>`__
* Batch system: `Slurm and HTCondor <https://submit.mit.edu/submit-users-guide/running.html>`__
* `Jupyter service <https://submit.mit.edu/jupyter/hub/spawn>`__
* `Filesystems <https://submit.mit.edu/submit-users-guide/storage.html>`__:

  * ``$HOME``: per-user directory, use only for inputs, source and scripts; backed up (10GB)
  * ``/work/submit/$USER``: community file system (100GB default)
  * ``/ceph/submit/data/user/<first_letter>/$USER``: storage for large files (1TB)

.. _building-submit-dependencies:

Dependencies
------------

Note that SubMIT does not provide software through system modules, as opposed to many other machines.
A recommended way to install packages is via Conda.


Install `Conda as per the documentation <https://submit.mit.edu/submit-users-guide/program.html#conda>`__ in your work area:

.. code-block:: bash

   cd /work/submit/$USER
   wget https://github.com/conda-forge/miniforge/releases/latest/download/Miniforge3-Linux-x86_64.sh
   bash Miniforge3-Linux-x86_64.sh


Follow the instructions on screen.
It is recommended to install Miniforge3 (and other software) in your work folder.
Therefore, when prompted ``Miniforge3 will now be installed into this location:``,
specify ``/work/submit/<username>/miniforge3`` for example.


eval "$(/home/submit/aforment/miniforge3/bin/conda shell.bash hook)"
conda init

Once Conda is installed adjust a few settings

.. code-block:: bash

   conda config --set auto_activate_base false
   conda config --add channels conda-forge
   conda config --set channel_priority strict

Then create a new WarpX environment with all the necessary dependencies:

.. code-block:: bash

   rm -rf $HOME/.conda/envs/warpx
   conda create --yes -n warpx -c conda-forge mamba conda-libmamba-solver
   conda activate warpx
   conda config --set solver libmamba
   mamba install --yes -c conda-forge python ipykernel ipympl h5py fast-histogram dask dask-jobqueue pyarrow blaspp boost ccache cmake compilers git lapackpp "openpmd-api=*=mpi_mpich*" openpmd-viewer packaging pytest python python-build make numpy pandas scipy setuptools yt "fftw=*=mpi_mpich*" pkg-config matplotlib mamba mpich mpi4py ninja pip virtualenv wheel libgomp


As the first step on future logins to SubMIT, remember to activate the ``warpx`` environment:

.. code-block:: bash

   conda activate warpx

Note that this is similar to what you could do on your :ref:``local machine <install-build-dependencies>``.


.. _building-submit-compilation:

Compilation
-----------

Use the following commands to download the WarpX source code:

.. code-block:: bash

   git clone https://github.com/BLAST-WarpX/warpx.git /work/submit/$USER/warpx


On SubMIT, you can run either on GPU nodes or CPU nodes.
There are three types of `available GPUs <https://submit.mit.edu/submit-users-guide/gpu.html>`__:

* ``nvidia_a30``: NVIDIA A30 GPUs
* ``Tesla_v100``: Tesla V100 GPUs
* ``nvidia_gtx1080``: NVIDIA GTX1080 GPUs


.. tab-set::

   .. tab-item:: GPUs Nodes

      To access one of the GPU partitions, use the following ``srun`` command or similar:

      .. code-block:: bash

         srun --partition=submit-gpu --gres=gpu:1 --cpus-per-gpu=8 --mem=64G --time=01:00:00 --pty bash


      To select, say, the NVIDIA A30 GPUs add ``--constraint=nvidia_a30`` to the ``srun`` command.

      Now that you have moved to the GPU partitions, the necessary NVIDIA software is readily available.
      Use the following :ref:`cmake commands <building-cmake>` to compile the application executable:

      .. code-block:: bash

         cd /work/submit/$USER/warpx
         rm -rf build_gpu_a30
         cmake -S . -B build_gpu_a30 -DWarpX_FFT=ON -DWarpX_DIMS=3 -DWarpX_COMPUTE=CUDA
         cmake --build build_gpu_a30 -j 8

      The WarpX application executables are now in ``/work/submit/$USER/warpx/build_gpu_a30/bin/``.
      Additionally, the following commands will install WarpX as a Python module:

      .. code-block:: bash

         cd /work/submit/$USER/warpx
         rm -rf build_gpu_a30_py

         cmake -S . -B build_gpu_a30_py -DWarpX_FFT=ON -DWarpX_DIMS=3 -DWarpX_COMPUTE=CUDA -DWarpX_APP=OFF -DWarpX_PYTHON=ON
         cmake --build build_gpu_a30_py -j 8 --target pip_install


   .. tab-item:: CPU Nodes

      Coming soon...



Now, you can :ref:`submit compute jobs <running-cpp-submit>` for WarpX :ref:`Python (PICMI) scripts <usage-picmi>` (:ref:`example scripts <usage-examples>`).
Or, you can use the WarpX executables to submit jobs (:ref:`example inputs <usage-examples>`).
For executables, you can reference their location in your :ref:`job script <running-cpp-submit>` or copy them to a location in simulation directory.


.. _running-cpp-submit:

Running
-------

.. tab-set::

   .. tab-item:: GPU Nodes

      The batch script below can be used to run a WarpX simulation on multiple nodes (change ``-N`` accordingly).
      This partition as up to `1536 nodes <https://docs.nersc.gov/systems/perlmutter/architecture/>`__.

      Replace descriptions between chevrons ``<>`` by relevant values, for instance ``<input file>`` could be ``plasma_mirror_inputs``.
      Note that we run one MPI rank per GPU.

      .. literalinclude:: ../../../../Tools/machines/submit-mit/submit_gpu.sbatch
         :language: bash
         :caption: You can copy this file from ``$HOME/src/warpx/Tools/machines/submit-mit/submit_gpu.sbatch``.

      To run a simulation, copy the lines above to a file ``submit_gpu.sbatch`` and run

      .. code-block:: bash

         sbatch perlmutter_gpu.sbatch

      to submit the job.


   .. tab-item:: CPU Nodes

      Coming soon...


.. _post-processing-submit:

Post-Processing
---------------

For post-processing, most users use `Jupyter service <https://submit.mit.edu/jupyter/hub/spawn>`__.

Create your own Jupyter kernel for post-processing:

.. code-block:: bash

   Coming soon


When opening a Jupyter notebook, just select ``WarpX`` from the list of available kernels on the top right of the notebook.

Additional software can be installed later on, e.g., in a Jupyter cell using ``!mamba install -y -c conda-forge ...``.
Software that is not available via conda can be installed via ``!python -m pip install ...``.
