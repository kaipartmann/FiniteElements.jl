# Newmark-Beta-Verfahren – Übersicht

Das Newmark-Beta-Verfahren ist ein numerisches Zeitintegrationsverfahren zur Lösung von Bewegungsgleichungen in der Strukturmechanik (insbesondere für dynamische Systeme wie z.B. Erdbebenanalysen).

## Ziel
- Zeitintegration von Bewegungsgleichungen der Form:  
  `M * ü(t) + C * ẋ(t) + K * x(t) = F(t)`
  - `M`: Massenmatrix  
  - `C`: Dämpfungsmatrix  
  - `K`: Steifigkeitsmatrix  
  - `x(t)`: Verschiebung  
  - `ẋ(t)`, `ü(t)`: Geschwindigkeit, Beschleunigung  
  - `F(t)`: äußere Kräfte  

## Grundprinzip
- Lösung erfolgt schrittweise über kleine Zeitschritte `Δt`
- Zustandsgrößen (`x`, `ẋ`, `ü`) werden für jeden Zeitschritt berechnet
- Zwei Parameter steuern das Verfahren:  
  - **β** (Beta): beeinflusst Stabilität und Genauigkeit  
  - **γ** (Gamma): steuert numerische Dämpfung

## Update-Formeln
Verschiebung und Geschwindigkeit bei `t + Δt` werden berechnet mit:
- `x_{n+1} = x_n + Δt * ẋ_n + Δt²/2 * ((1 - 2β) * ü_n + 2β * ü_{n+1})`
- `ẋ_{n+1} = ẋ_n + Δt * ((1 - γ) * ü_n + γ * ü_{n+1})`

## Typische Werte für die Parameter
| Verfahren         | β       | γ       | Eigenschaften               |
|------------------|---------|---------|-----------------------------|
| Explizit         | 0       | 0.5     | Nicht stabil für alle Δt    |
| Konservativ      | 1/6     | 0.5     | Energieerhaltend            |
| Standard (impl.) | 1/4     | 0.5     | Unkond. stabil, kein Drift  |

## Vorteile
- Einfache Implementierung
- Flexibel durch Wahl von β und γ
- Breite Anwendung in FEM und Strukturmechanik

## Nachteile
- Kann numerische Dämpfung verursachen
- Bei falscher Wahl von β/γ instabil oder ungenau

## Anwendung
- Zeitbereichsanalyse in der Finite-Elemente-Methode (FEM)
- Simulation von Erdbebenlasten
- Schwingungsprobleme in Maschinenbau, Bauwesen etc.
