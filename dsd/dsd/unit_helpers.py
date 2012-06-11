import numpy as np

def no_units(num):
    try:
        return num.simplified.magnitude
    except AttributeError:
        return num

try:
    import functools
    import quantities as pq
    from quantities.decorators import quantitizer as quant

    def quantitizer(handler):
        def dec(func):
            return quant(func, handler)
        return dec

    # Make a decorator to try to do appropriate unit conversion
    def check_units(**units):
        def decorator(func):
            @functools.wraps(func)
            def wrapper(*args, **kwargs):
                args = list(args)

                # By default, we do not want units back out
                returnUnits = False

                # Handle keyword args
                for kw,val in kwargs.items():
                    if isinstance(val, pq.Quantity):
                        # Since one was passed in, return a Quantity
                        returnUnits = True
                    else:
                        if kw in units:
                            kwargs[kw] = pq.Quantity(val, units=units[kw])

                # Now handle positional args by linking them in the order
                # given to the names in the code object.
                for argInd,(arg,val) in enumerate(
                        zip(func.func_code.co_varnames, args)):
                    if isinstance(val, pq.Quantity):
                        returnUnits = True
                    else:
                        if arg in units:
                            args[argInd] = pq.Quantity(val, units=units[arg])

                # Call the function
                ret = func(*args, **kwargs)

                # Make sure the returned type is appropriate
                if not returnUnits:
                    ret = ret.simplified.magnitude
                return ret
            return wrapper
        return decorator

    # Make a decorator to force the proper units if we're given a
    # quantity
    def force_units(return_units, **units):
        def decorator(func):
            @functools.wraps(func)
            def wrapper(*args, **kwargs):
                args = list(args)

                # By default, we do not want units back out
                returnUnits = False

                # Handle keyword args
                for kw,val in kwargs.items():
                    if isinstance(val, pq.Quantity):
                        # Since one was passed in, return a Quantity
                        returnUnits = True
                        newUnits = units.get(kw, '')
                        kwargs[kw] = val.rescale(newUnits).magnitude

                # Now handle positional args by linking them in the order
                # given to the names in the code object.
                for argInd,(arg,val) in enumerate(
                        zip(func.func_code.co_varnames, args)):
                    if isinstance(val, pq.Quantity):
                        returnUnits = True
                        newUnits = units.get(arg, '')
                        args[argInd] = val.rescale(newUnits).magnitude

                # Call the function
                ret = func(*args, **kwargs)

                # Force return type
                if returnUnits:
                    try:
                        if return_units is not None:
                            ret = pq.Quantity(ret, units=return_units)
                    except TypeError:
                        ret = tuple(pq.Quantity(r, units=u) for u,r
                                in zip(return_units, ret))
                return ret
            return wrapper
        return decorator

    # Needed to we can call exp with unitless quantities that just need
    # proper scaling first, like m/cm
    def exp(val, *args, **kwargs):
        try:
            val = val.simplified
        except AttributeError:
            pass
        return np.exp(val, *args, **kwargs)

    def angle(val, deg=False):
        try:
            val = val.simplified
            ret = np.angle(val, deg)
            if deg:
                return pq.Quantity(ret, units='degrees')
            else:
                return pq.Quantity(ret, units='radians')
        except AttributeError:
            return np.angle(val, deg)

    unit_dict = dict()
    unit_dict['density_water'] = pq.kilogram / pq.meter**3
    unit_dict['mp_N0'] = pq.meter**-4
    def update_consts(local_dict):
        for kw in unit_dict:
            if kw in local_dict:
                local_dict[kw] = local_dict[kw] * unit_dict[kw]

except ImportError:
    def check_units(**kwargs):
        def dec(func):
            return func
        return dec

    def quantitizer(handler):
        def dec(func):
            return func
        return dec

    force_units = check_units

    exp = np.exp
    angle = np.angle

    def update_consts(local_dict):
        pass
