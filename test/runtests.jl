# Copyright 2020 Digital Domain 3.0
#
# Licensed under the Apache License, Version 2.0 (the "Apache License")
# with the following modification; you may not use this file except in
# compliance with the Apache License and the following modification to it:
# Section 6. Trademarks. is deleted and replaced with:
#
# 6. Trademarks. This License does not grant permission to use the trade
#    names, trademarks, service marks, or product names of the Licensor
#    and its affiliates, except as required to comply with Section 4(c) of
#    the License and to reproduce the content of the NOTICE file.
#
# You may obtain a copy of the Apache License at
#
#     http://www.apache.org/licenses/LICENSE-2.0
#
# Unless required by applicable law or agreed to in writing, software
# distributed under the Apache License with the above modification is
# distributed on an "AS IS" BASIS, WITHOUT WARRANTIES OR CONDITIONS OF ANY
# KIND, either express or implied. See the Apache License for the specific
# language governing permissions and limitations under the Apache License.

using Test
using Blades
using KVectors
using LinearAlgebra

module PG3
  using Blades
  @generate_basis("0+++")
end

module G2
  using Blades
  # need last param to be true for differential forms
  @generate_basis("++",false,true,true)
end

module G3
  using Blades
  @generate_basis("+++",false,true,true)
end

module G5
  using Blades
  @generate_basis("+++++",false,true,true)
end

using .PG3
using .G2
using .G3
using .G5
  
@testset "KVectors" begin
  e₁, e₂, e₃, e₄ = alle( PG3, 4)[1:4]
  a = e₁(1.0); b = e₂(2.0); c = e₃(3.0); d = e₄(4.0)
  @test typeof(a+a) <: Blade
  @test (typeof(a+b) <: Blade) == false
  @test typeof(a+b) <: KVector
  B = a+b
  B2 = b+c
  @test a+b+c+b+a == e₁(2.0) + e₂(4.0) + e₃(3.0) == 2.0a+2.0b+c
  @test B+B2 == B2+B == a+2.0b+c
  @test (a+b)+(c+d) == a+b+c+d
  @test 5.0*(a+b+c+d) == 5.0*b+5.0*c+5.0*d+5.0*a
  @test grade(a+b+c) == 1
  @test grade(KVector(a*b)) == 2
  @test dual(dual(B)) == B
  @test zero(KVector(a))+b == KVector(b) == b + zero(KVector(a))
  @test zero(B)*B == zero(B) == B*zero(B) == B*zero(KVector{Float64, grade(B), 1})
  @test 3.0*B == B+B+B == B*3.0

end

