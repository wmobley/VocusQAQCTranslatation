pro test_compound_class

  compile_opt idl2
  @unittest_setup

  ; Test case for organic compound with carbon and oxygen
  result = CompoundClass('C2H4O2')
  assert, array_equal(result, ['org_', 'HCO_', 'HCO2', 'unsat=1.0C', 'OSC=-1.00000', 'nC= 2', 'nH= 4', 'nO= 2', 'nN= 0']), 'Failed for C2H4O2'

  ; Test case for organic compound with nitrogen
  result = CompoundClass('C3H7NO2')
  assert, array_equal(result, ['org_', 'orgN_', 'orgN1', 'NHCO2', 'unsat=1.0N', 'OSC=-1.33333', 'nC= 3', 'nH= 7', 'nO= 2', 'nN= 1']), 'Failed for C3H7NO2'

  ; Test case for inorganic compound
  result = CompoundClass('H2SO4')
  assert, array_equal(result, ['inorg', 'unsat=0.0C', 'OSC=-9999.00', 'nC= 0', 'nH= 2', 'nO= 4', 'nN= 0']), 'Failed for H2SO4'

  ; Test case for hydrocarbon
  result = CompoundClass('C6H14')
  assert, array_equal(result, ['org_', 'HC_', 'unsat=0.0C', 'OSC=-2.33333', 'nC= 6', 'nH=14', 'nO= 0', 'nN= 0']), 'Failed for C6H14'

  ; Test case for compound with 13C
  result = CompoundClass('13CC5H12')
  assert, array_equal(result, ['org_', 'HC_', 'unsat=0.0C', 'OSC=-2.33333', 'nC= 6', 'nH=12', 'nO= 0', 'nN= 0']), 'Failed for 13CC5H12'

  ; Test case for organic sulfate
  result = CompoundClass('C2H5SO4')
  assert, array_equal(result, ['org_', 'orgSO4', 'unsat=0.5C', 'OSC=1.00000', 'nC= 2', 'nH= 5', 'nO= 4', 'nN= 0']), 'Failed for C2H5SO4'

  @unittest_teardown

end
