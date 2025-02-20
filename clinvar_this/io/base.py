"""Base module supporting for I/O of various formats to define submissions."""

from abc import ABC, abstractmethod
import pathlib
import typing

from pydantic import BaseModel, ConfigDict

from clinvar_api.models import SubmissionContainer
from clinvar_api.msg.sub_payload import AlleleOrigin, CollectionMethod
from clinvar_this import exceptions


class BatchMetadata(BaseModel):
    """Batch-wide settings for import.

    The properties will be assigned to all variants/samples in the batch.
    """

    model_config = ConfigDict(frozen=True)

    collection_method: typing.Optional[CollectionMethod] = None
    allele_origin: typing.Optional[AlleleOrigin] = None


class TransformIO(ABC):
    """Base class for transforming input data from various formats into submission format"""

    batch_metadata_defaults: typing.Dict[str, typing.Any]
    batch_metadata_model: BatchMetadata

    def batch_metadata_from_mapping(
        self,
        keys_values: typing.Iterable[str],
        *,
        use_defaults: bool = False,
    ) -> BatchMetadata:
        """Convert configuration from ``KEY=VALUE`` strings to ``batch_metadata_model``

        Default values can be used (should be on import but not on update).

        :param keys_values: Key-values for ``batch_metadata_model``
        :param use_defaults: Whether to use the key-value pairs from
            ``batch_metadata_defaults``
        :return: Batch metadata model
        """
        kwargs = {}
        if use_defaults:
            for key, value in self.batch_metadata_defaults.model_dump().items():
                kwargs.setdefault(key, value)
        else:
            field_types = {
                name: value
                for (name, value) in typing.get_type_hints(
                    self.batch_metadata_model
                ).items()
            }
            for key_value in keys_values:
                if "=" not in key_value:
                    raise exceptions.ArgumentsError(
                        f"Invalid key/value pair in {key_value}"
                    )
                key, value = key_value.split("=")
                if key in field_types:
                    try:
                        # We need to ignore types as mypy 1.6.0 warns for "expected type[any] but
                        # found "type[any] | None".
                        kwargs[key] = field_types[key].model_validate(value)  # type: ignore
                    except ValueError:
                        raise exceptions.ArgumentsError(
                            f"Failed to parse {value} as for key {key}"
                        )

        return self.batch_metadata_model(**kwargs)

    @abstractmethod
    def _read_file(inputf: typing.TextIO) -> typing.List:
        """Read input from a given file and transform to required input data structure

        :param inputf: Text file-like object containing input data
        :return: A list of structured records
        """

    @abstractmethod
    def records_to_submission_container(
        self, *args, **kwargs
    ) -> typing.List[SubmissionContainer]:
        """Transform structured records to submission container data structures

        :return: A list of submission container data structures
        """

    def read_file(
        self,
        file: typing.Optional[typing.TextIO] = None,
        path: typing.Union[None, str, pathlib.Path] = None,
    ) -> typing.List[dict]:
        """Read input from either text file-like object or file path

        :param file: Text file-like object containing input data
        :param path: The path to the input file
        :return: A list of dictionaries representing the transformed input data
        """
        if file:
            return self._read_file(file)
        elif path:
            with pathlib.Path(path).open("rt") as inputf:
                return self._read_file(inputf)
        else:
            raise TypeError("You have to provide either file or path")
