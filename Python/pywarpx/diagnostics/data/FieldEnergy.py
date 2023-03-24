from BaseReader import DataReader
import numpy as np
import pandas as pd


class FieldEnergyData(DataReader):
    """
    Reader for the FieldEnergy reduced diagnostic
    """

    def __init__(self, run_directory, file_prefix='FieldEnergy'):
        """
        Parameters
        ----------
        run_directory: string
            Path to the run directory of WarpX.
        file_prefix : string
            Name of the file containing the current reduced diagnostic data.
        """

        super().__init__(run_directory)

        self.data_file_prefix = file_prefix
        self.data_file_suffix = '.txt'

        self.has_species = False

    def get_nlev(self):
        """Returns the number of mesh refinement levels in the simulation."""

        return int((self.get_ncols() - 2.)/3.)

    def get_valid_args(self):
        nlev = self.get_nlev()
        valid_args = [s for n in range(nlev) for s in ['total_lev'+str(n), 'E_lev'+str(n), 'B_lev'+str(n)]]
        return valid_args


    def get_data(self, *args, steps=None, times=None):
        """
        Arguments:
            string:
                Can be one or more among the valid arguments.
        Keyword arguments:
            steps = list or np.array of integers or None (optional):
                Timesteps at which the desidered output will be returned.
                If equal to None or not specified then all timesteps are given.
            times = list or np.array of numbers or None (optional):
                The desidered output will be returned at the closest availble times.
                If equal to None or not specified then all timesteps are given
        Output:
            pandas dataframe with columns: steps, times, requested data
        """

        valid_args = self.get_valid_args()

        # get data using parent class
        data = DataReader.get_data(self, valid_args, steps, times, *args)

        # restrict to the desidered steps or times
        restricted_data = self.restrict_data(data, steps, times)

        # rename columns
        restricted_data.columns = ['steps', 'times', *[name for name in args]]

        return restricted_data
        
        
