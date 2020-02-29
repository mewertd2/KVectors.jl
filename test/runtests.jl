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
using StaticArrays

module PG3
  using Blades
  @generate_basis("0+++")
end
using .PG3

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
  @test iszero(B∧B)
  @test 2.0∧B∧B2 == B∧b∧2.0 - c∧B∧2.0
  @test grade(B∧B2) == 2
  @test grade(B,1) == B
  @test reverse(B) == B
  @test reverse(B∧B2) == -B∧B2
  @test B-B == -a - -B - b

  @test dual(prune(B-B)) == prune(B-B)
  @test iszero(KVector(0.0e₂))
  @test iszero(1.0e₁-1.0e₁)
  @test dual(B) == !B

#!me passes with Blades v0.1.1+  @test normalize(KVector(-2.2(e₂∧e₃))) == KVector(normalize(-2.2(e₂∧e₃))) == KVectors.normalize_safe(KVector(-2.2(e₂∧e₃)))
end

module G3
  using Blades
  @generate_basis("+++",false,true,true)
end
using .G3
 
@testset "More KVectors" begin
  e₁, e₂, e₃ = alle(G3, 3)[1:3]
  𝐼 = alle(G3,3)[end]

  a = sortbasis(1.0e₁ + 3.0e₃)

  @test a == KVector([1.0,0.0,3.0]) |> KVectors.prune
  @test first(a) == a[1]
  @test firstindex(a) == 1
  @test a[end] == a[2]
  @test isnull(a) == false
  @test length(a) == 2
  @test isempty(a) == false
  @test [i for i in a] == map(i->i, a) == (i for i in a) |> collect 
  B = -1.0(e₁∧e₂) + 2.0(e₁∧e₃)
  x = 0.0
  for i in a
    x = x+scalar(i)
  end
  @test x == mapreduce(scalar, +, a)

  @test conj(a) == a
  @test conj(B) == -B
  @test KVector(a) == a
  @test KVector([1,2,3], 𝐼) == 1e₁+2e₂+3e₃
  @test pseudoscalar(a) == pseudoscalar(a[1])
  @test grade(⟂(a)∧a) == grade(pseudoscalar(a))
  @test a/2.0 == a*0.5

  @test coords(a) == scalar.(sortbasis(a+0.0e₂))
  @test coords(a[1]) == [scalar(a[1]), 0.0, 0.0]
  @test KVectors.prune(KVector(coords(a) .* basis_1blades(a))) == a
  @test norm(basis_1vector(a)) == sqrt(3.0)
  @test KVectors.norm_sqr(a) == mapreduce(aᵢ->aᵢ*aᵢ, +, a)
  @test norm(KVectors.normalize_safe(a)) == norm(normalize(a))

  # warning this relys on LinearAlebra:⋅ 
  # which only works on KVectors of grade 1 with matching sorted basis 1-vectors present
  # for a proper inner product you should be using Multivectors.jl
  @test det(1.0e₁+1.0e₂, 1.0e₁+1.0e₂) == 1.0
end
