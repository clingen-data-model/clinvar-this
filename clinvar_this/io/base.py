"""Base module supporting for I/O of various formats to define submissions."""

from abc import ABC, abstractmethod
import pathlib
import typing


from clinvar_api.models import SubmissionContainer


class TransformIO(ABC):
    """Base class for transforming input data from various formats into submission format"""

    @abstractmethod
    def _read_file(self, inputf: typing.TextIO) -> typing.List:
        """Read input from a given file and transform to required input data structure

        :param inputf: Text file-like object containing input data
        :return: A list of structured records
        """

    @abstractmethod
    def records_to_submission_container(
        self, *args, **kwargs
    ) -> SubmissionContainer:
        """Transform structured records to submission container data structures

        :return: Submission container data structures
        """

    def read_file(
        self,
        file: typing.Optional[typing.TextIO] = None,
        path: typing.Union[None, str, pathlib.Path] = None,
    ) -> typing.List:
        """Read input from either text file-like object or file path

        :param file: Text file-like object containing input data
        :param path: The path to the input file
        :return: A list of the transformed input data
        """
        if file:
            return self._read_file(file)
        elif path:
            with pathlib.Path(path).open("rt") as inputf:
                return self._read_file(inputf)
        else:
            raise TypeError("You have to provide either file or path")
