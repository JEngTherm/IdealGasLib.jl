#----------------------------------------------------------------------------------------------#
#                                         interface.jl                                         #
#----------------------------------------------------------------------------------------------#

# Interface (bare) functions

"""
Function to return a particular substance's name.
"""
function name end

"""
Function to return a particular substance's chemical formula.
"""
function form end

"""
Function to return a particular substance's reference state's temperature.
"""
function Tref end

"""
Function to return a particular substance's reference state's specific entropy.
"""
function sref end


# Interface exports
export name, form, Tref, sref


