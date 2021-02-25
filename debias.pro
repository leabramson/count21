function debias, structure

  norms = scaleToMean(structure)
  structure.COUNTS *= norms.GLOBAL

  return, structure
end
