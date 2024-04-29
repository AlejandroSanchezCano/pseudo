# Built-in modules
from typing import Any

# Third-party modules
import pickle

def pickling(data: Any, path: str) -> None:
    '''
    Pickle an object and store it.

    Parameters
    ----------
    data : Any
        Pickable object that will be stored.
    path : str
        Storing path.
    '''
    with open(path, 'wb') as handle:
        pickle.dump(
            obj = data,
            file = handle, 
            protocol = pickle.HIGHEST_PROTOCOL
            )

def unpickling(path: str) -> Any:
    '''
    Retrieves and unpickles a pickled object.

    Parameters
    ----------
    path : str
        Storing path of the object to unpickle.

    Returns
    -------
    Any
        Unpickled object.
    '''
    with open(path, 'rb') as handle:
        return pickle.load(file = handle)