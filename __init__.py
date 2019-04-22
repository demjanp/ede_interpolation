# -*- coding: utf-8 -*-
"""
/***************************************************************************
 EDEInterpolation
                                 A QGIS plugin
 Evidence Density Estimation (EDE) interpolation of archaeological settlement data.
 (Demján and Dreslerová 2016) https://doi.org/10.1016/j.jas.2016.04.003
                             -------------------
        begin                : 2019-04-10
        copyright            : (C) 2019 by Peter Demján
        email                : peter.demjan@gmail.com
        git sha              : $Format:%H$
 ***************************************************************************/

/***************************************************************************
 *                                                                         *
 *   This program is free software; you can redistribute it and/or modify  *
 *   it under the terms of the GNU General Public License as published by  *
 *   the Free Software Foundation; either version 2 of the License, or     *
 *   (at your option) any later version.                                   *
 *                                                                         *
 ***************************************************************************/
 This script initializes the plugin, making it known to QGIS.
"""


# noinspection PyPep8Naming
def classFactory(iface):  # pylint: disable=invalid-name
    """Load EDEInterpolation class from file EDEInterpolation.

    :param iface: A QGIS interface instance.
    :type iface: QgsInterface
    """
    #
    from .ede_interpolation import EDEInterpolation
    return EDEInterpolation(iface)
