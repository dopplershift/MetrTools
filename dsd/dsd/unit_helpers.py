import numpy as np

try:
    import functools
    import quantities as pq

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
                        newUnits = units.get(kw, '')
                        kwargs[kw] = pq.Quantity(val, units=units[kw])

                # Now handle positional args by linking them in the order
                # given to the names in the code object.
                for argInd,(arg,val) in enumerate(
                        zip(func.func_code.co_varnames, args)):
                    if isinstance(val, pq.Quantity):
                        returnUnits = True
                    else:
                        newUnits = units.get(arg, '')
                        args[argInd] = pq.Quantity(val, units=newUnits)

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
                    ret = pq.Quantity(ret, units=return_units)
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

    unit_dict = dict()
    unit_dict['density_water'] = pq.kilogram / pq.meter**3
    unit_dict['mp_N0'] = pq.meter**-4
    unit_dict['gamma_power_scale'] = pq.meter
    def update_consts(local_dict):
        for kw in unit_dict:
            if kw in local_dict:
                local_dict[kw] = local_dict[kw] * unit_dict[kw]

except ImportError:
    def check_units(**kwargs):
        def dec(func):
            return func
        return dec

    force_units = check_units

    exp = np.exp

    def update_consts(local_dict):
        pass
