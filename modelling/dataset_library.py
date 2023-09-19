import myokit


class ProtocolLibrary(object):
    """
    A library class with known protocols
    """
    def __init__(self):
        super(ProtocolLibrary, self).__init__()

    def Milnes(self, t_max):
        protocol = myokit.Protocol()
        protocol.schedule(-80, 0, 800, period=t_max)
        protocol.schedule(-90, 800, 100, period=t_max)
        protocol.schedule(-80, 900, 100, period=t_max)
        protocol.schedule(-80, 11000, 14000, period=t_max)

        return protocol

    def current_impulse(self, t_max, offset=50):
        return myokit.pacing.blocktrain(t_max, 1, offset=offset)

    def validation3(self, t_max):
        protocol = myokit.Protocol()
        protocol.schedule(-80, 0, 100)  # , period=t_max)
        protocol.schedule(-40, 100, 50)  # , period=t_max)
        protocol.schedule(20, 150, 500)  # , period=t_max)
        protocol.schedule(-40, 650, 500)  # , period=t_max)
        protocol.schedule(-80, 1150, 200)  # , period=t_max)

        return protocol

    def hERG_validation(self, t_max):
        protocol = myokit.Protocol()
        protocol.schedule(-80, 0, 100)
        protocol.schedule(20, 100, 900)
        protocol.schedule(-80, 1000, 1000)

        return protocol

    def Pneg80(self, t_max):
        protocol = myokit.Protocol()
        protocol.schedule(-80, 0, 200, period=t_max)
        protocol.schedule(20, 200, 500, period=t_max)
        protocol.schedule(-50, 700, 200, period=t_max)
        protocol.schedule(-80, 900, 4500, period=t_max)

        return protocol

    def P0(self, t_max):
        protocol = myokit.Protocol()
        protocol.schedule(-80, 0, 100, period=t_max)
        protocol.schedule(-60, 5100, 200, period=t_max)
        protocol.schedule(-80, 5300, 100, period=t_max)

        return protocol

    def P40(self, t_max):
        protocol = myokit.Protocol()
        protocol.schedule(-80, 0, 100, period=t_max)
        protocol.schedule(40, 100, 5000, period=t_max)
        protocol.schedule(-60, 5100, 200, period=t_max)
        protocol.schedule(-80, 5300, 100, period=t_max)

        return protocol
